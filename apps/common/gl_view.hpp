#pragma once

#include <GL/glew.h>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "font_data.hpp"

namespace spp::app {
using Color=std::array<float,4>;
inline constexpr Color ink{0.86f,0.88f,0.91f,1}, muted{0.48f,0.53f,0.60f,1}, accent{0.90f,0.72f,0.40f,1};
struct Rect {
    float x,y,w,h;
    bool contains(double px,double py) const {
        return px>=static_cast<double>(x) && px<static_cast<double>(x+w)
            && py>=static_cast<double>(y) && py<static_cast<double>(y+h);
    }
};

inline GLuint shader(GLenum type,const char* source) {
    const GLuint id=glCreateShader(type);
    glShaderSource(id,1,&source,nullptr); glCompileShader(id);
    GLint good=0; glGetShaderiv(id,GL_COMPILE_STATUS,&good);
    if (!good) {
        std::array<char,4096> log{}; glGetShaderInfoLog(id,static_cast<GLsizei>(log.size()),nullptr,log.data());
        glDeleteShader(id); throw std::runtime_error(std::string("OpenGL shader: ")+log.data());
    }
    return id;
}
inline GLuint program(const char* vertex,const char* fragment) {
    const GLuint v=shader(GL_VERTEX_SHADER,vertex);
    GLuint f=0;
    try { f=shader(GL_FRAGMENT_SHADER,fragment); } catch (...) { glDeleteShader(v); throw; }
    const GLuint p=glCreateProgram(); glAttachShader(p,v); glAttachShader(p,f); glLinkProgram(p);
    glDeleteShader(v); glDeleteShader(f);
    GLint good=0; glGetProgramiv(p,GL_LINK_STATUS,&good);
    if (!good) {
        std::array<char,4096> log{}; glGetProgramInfoLog(p,static_cast<GLsizei>(log.size()),nullptr,log.data());
        glDeleteProgram(p); throw std::runtime_error(std::string("OpenGL link: ")+log.data());
    }
    return p;
}

class View {
    struct Vertex { float x,y,u,v; Color color; };
    GLuint ui_=0,disk_=0,vao_=0,vbo_=0,font_=0;
    GLint ui_viewport_=0,disk_viewport_=0,disk_center_=0,disk_radius_=0,light_=0,earthshine_=0,area_=0;
    std::vector<Vertex> vertices_;
    float width_=1280,height_=900;
public:
    View() {
        constexpr auto uv=R"(#version 330 core
layout(location=0) in vec2 position; layout(location=1) in vec2 tex;
layout(location=2) in vec4 color; uniform vec2 viewport;
out vec2 uv; out vec4 tint;
void main(){gl_Position=vec4(position.x/viewport.x*2-1,1-position.y/viewport.y*2,0,1);uv=tex;tint=color;}
)";
        constexpr auto uf=R"(#version 330 core
in vec2 uv; in vec4 tint; uniform sampler2D atlas; out vec4 out_color;
void main(){out_color=vec4(tint.rgb,tint.a*(uv.x<0?1:texture(atlas,uv).r));}
)";
        constexpr auto dv=R"(#version 330 core
void main(){vec2 p=vec2((gl_VertexID<<1)&2,gl_VertexID&2);gl_Position=vec4(p*2-1,0,1);}
)";
        constexpr auto df=R"(#version 330 core
uniform vec2 viewport,center; uniform float radius,earthshine; uniform vec3 light;
uniform bool area; out vec4 out_color;
void main(){
    vec2 p=(vec2(gl_FragCoord.x,viewport.y-gl_FragCoord.y)-center)/radius;
    p.y=-p.y; float r2=dot(p,p); if(r2>1)discard;
    vec3 n=vec3(p,sqrt(max(0,1-r2))); float mu0=max(dot(n,light),0);
    float lit=area?(mu0>0?0.88:0):0.26*mu0+0.74*mu0/max(mu0+n.z,0.000001);
    lit=max(lit,earthshine*(0.6+0.4*n.z));
    vec3 albedo=vec3(0.84,0.85,0.87); vec3 c=albedo*pow(max(lit,0),0.454545);
    float coverage=1-smoothstep(1-1.5/radius,1,sqrt(r2)); out_color=vec4(c,coverage);
}
)";
        ui_=program(uv,uf); disk_=program(dv,df);
        ui_viewport_=glGetUniformLocation(ui_,"viewport");
        disk_viewport_=glGetUniformLocation(disk_,"viewport"); disk_center_=glGetUniformLocation(disk_,"center");
        disk_radius_=glGetUniformLocation(disk_,"radius"); light_=glGetUniformLocation(disk_,"light");
        earthshine_=glGetUniformLocation(disk_,"earthshine"); area_=glGetUniformLocation(disk_,"area");
        glGenVertexArrays(1,&vao_); glGenBuffers(1,&vbo_);
        glBindVertexArray(vao_); glBindBuffer(GL_ARRAY_BUFFER,vbo_);
        glEnableVertexAttribArray(0); glVertexAttribPointer(0,2,GL_FLOAT,GL_FALSE,sizeof(Vertex),nullptr);
        glEnableVertexAttribArray(1); glVertexAttribPointer(1,2,GL_FLOAT,GL_FALSE,sizeof(Vertex),reinterpret_cast<void*>(2*sizeof(float)));
        glEnableVertexAttribArray(2); glVertexAttribPointer(2,4,GL_FLOAT,GL_FALSE,sizeof(Vertex),reinterpret_cast<void*>(4*sizeof(float)));
        std::vector<unsigned char> pixels; pixels.reserve(font::width*font::height);
        for (std::size_t i=0;i<sizeof(font::rle);i+=2) pixels.insert(pixels.end(),font::rle[i],font::rle[i+1]);
        if (pixels.size()!=font::width*font::height) throw std::runtime_error("invalid embedded font atlas");
        glGenTextures(1,&font_); glBindTexture(GL_TEXTURE_2D,font_); glPixelStorei(GL_UNPACK_ALIGNMENT,1);
        glTexImage2D(GL_TEXTURE_2D,0,GL_R8,font::width,font::height,0,GL_RED,GL_UNSIGNED_BYTE,pixels.data());
        glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR); glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_CLAMP_TO_EDGE); glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_CLAMP_TO_EDGE);
        vertices_.reserve(24000);
    }
    View(const View&)=delete; View& operator=(const View&)=delete;
    ~View() { glDeleteProgram(ui_);glDeleteProgram(disk_);glDeleteBuffers(1,&vbo_);glDeleteVertexArrays(1,&vao_);glDeleteTextures(1,&font_); }
    void begin(int w,int h) {
        width_=static_cast<float>(w);height_=static_cast<float>(h);vertices_.clear();
        glViewport(0,0,w,h);glDisable(GL_DEPTH_TEST);glEnable(GL_BLEND);glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
        glClearColor(0.035f,0.043f,0.057f,1);glClear(GL_COLOR_BUFFER_BIT);glBindVertexArray(vao_);
    }
    void quad(Rect r,Color color,float u0=-1,float v0=-1,float u1=-1,float v1=-1) {
        const Vertex a{r.x,r.y,u0,v0,color}, b{r.x+r.w,r.y,u1,v0,color};
        const Vertex c{r.x+r.w,r.y+r.h,u1,v1,color}, d{r.x,r.y+r.h,u0,v1,color};
        vertices_.insert(vertices_.end(),{a,b,c,a,c,d});
    }
    void line(float ax,float ay,float bx,float by,Color color,float thickness=1) {
        const float len=std::hypot(bx-ax,by-ay); if (len<=0)return;
        const float dx=(by-ay)*0.5f*thickness/len,dy=(ax-bx)*0.5f*thickness/len;
        const Vertex a{ax-dx,ay-dy,-1,-1,color}, b{ax+dx,ay+dy,-1,-1,color};
        const Vertex c{bx+dx,by+dy,-1,-1,color}, d{bx-dx,by-dy,-1,-1,color};
        vertices_.insert(vertices_.end(),{a,b,c,a,c,d});
    }
    void text(float x,float y,std::string_view s,float scale=0.85f,Color color=ink) {
        const float start=x;
        for (unsigned char c:s) {
            if (c=='\n') {x=start;y+=font::cell_height*scale;continue;}
            const auto i=static_cast<unsigned>(c>=32 && c<127?c-32:'?'-32);
            const float u0=static_cast<float>((i%16)*font::cell_width)/font::width;
            const float v0=static_cast<float>((i/16)*font::cell_height)/font::height;
            quad({x,y,font::cell_width*scale,font::cell_height*scale},color,u0,v0,
                u0+static_cast<float>(font::cell_width)/font::width,v0+static_cast<float>(font::cell_height)/font::height);
            x+=font::advance[i]*scale;
        }
    }
    void disk(float x,float y,float radius,const std::array<float,3>& light,bool earthshine=false,bool area=false) {
        glUseProgram(disk_);glUniform2f(disk_viewport_,width_,height_);glUniform2f(disk_center_,x,y);
        glUniform1f(disk_radius_,radius);glUniform3fv(light_,1,light.data());
        glUniform1f(earthshine_,earthshine?0.012f:0.0f);glUniform1i(area_,area);
        // Clip the full-screen triangle to the disk's bounding rectangle.
        const int left=static_cast<int>(std::floor(x-radius-2)),bottom=static_cast<int>(std::floor(height_-y-radius-2));
        glEnable(GL_SCISSOR_TEST);glScissor(left,bottom,static_cast<int>(2*radius+4),static_cast<int>(2*radius+4));
        glDrawArrays(GL_TRIANGLES,0,3);glDisable(GL_SCISSOR_TEST);
    }
    void end() {
        glUseProgram(ui_);glUniform2f(ui_viewport_,width_,height_);glActiveTexture(GL_TEXTURE0);glBindTexture(GL_TEXTURE_2D,font_);
        glBindBuffer(GL_ARRAY_BUFFER,vbo_);glBufferData(GL_ARRAY_BUFFER,static_cast<GLsizeiptr>(vertices_.size()*sizeof(Vertex)),vertices_.data(),GL_STREAM_DRAW);
        glDrawArrays(GL_TRIANGLES,0,static_cast<GLsizei>(vertices_.size()));
    }
};
}  // namespace spp::app
