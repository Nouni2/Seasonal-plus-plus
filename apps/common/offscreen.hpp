#pragma once

#include <stdexcept>
#if SPP_MOON_EGL
#include <EGL/egl.h>
#include <EGL/eglext.h>
#endif

namespace spp::app {
// Used only for explicit rendering performance runs on a headless machine.
class Offscreen {
#if SPP_MOON_EGL
    EGLDisplay display_=EGL_NO_DISPLAY;
    EGLSurface surface_=EGL_NO_SURFACE;
    EGLContext context_=EGL_NO_CONTEXT;
    void release() {
        if (display_!=EGL_NO_DISPLAY) {
            eglMakeCurrent(display_,EGL_NO_SURFACE,EGL_NO_SURFACE,EGL_NO_CONTEXT);
            if (context_!=EGL_NO_CONTEXT) eglDestroyContext(display_,context_);
            if (surface_!=EGL_NO_SURFACE) eglDestroySurface(display_,surface_);
            eglTerminate(display_);display_=EGL_NO_DISPLAY;
        }
    }
#endif
public:
    Offscreen(int w,int h) {
#if SPP_MOON_EGL
        try {
            auto get_platform=reinterpret_cast<PFNEGLGETPLATFORMDISPLAYEXTPROC>(eglGetProcAddress("eglGetPlatformDisplayEXT"));
            if (!get_platform) throw std::runtime_error("EGL platform display extension unavailable");
            display_=get_platform(EGL_PLATFORM_SURFACELESS_MESA,EGL_DEFAULT_DISPLAY,nullptr);
            if (display_==EGL_NO_DISPLAY || !eglInitialize(display_,nullptr,nullptr) || !eglBindAPI(EGL_OPENGL_API))
                throw std::runtime_error("cannot initialize surfaceless EGL OpenGL");
            const EGLint attributes[]={EGL_SURFACE_TYPE,EGL_PBUFFER_BIT,EGL_RENDERABLE_TYPE,EGL_OPENGL_BIT,
                EGL_RED_SIZE,8,EGL_GREEN_SIZE,8,EGL_BLUE_SIZE,8,EGL_ALPHA_SIZE,8,EGL_NONE};
            EGLConfig config{};EGLint count=0;
            if (!eglChooseConfig(display_,attributes,&config,1,&count) || !count)
                throw std::runtime_error("no EGL OpenGL pbuffer configuration");
            const EGLint size[]={EGL_WIDTH,w,EGL_HEIGHT,h,EGL_NONE};
            surface_=eglCreatePbufferSurface(display_,config,size);
            const EGLint version[]={EGL_CONTEXT_MAJOR_VERSION,3,EGL_CONTEXT_MINOR_VERSION,3,
                EGL_CONTEXT_OPENGL_PROFILE_MASK,EGL_CONTEXT_OPENGL_CORE_PROFILE_BIT,EGL_NONE};
            context_=eglCreateContext(display_,config,EGL_NO_CONTEXT,version);
            if (surface_==EGL_NO_SURFACE || context_==EGL_NO_CONTEXT || !eglMakeCurrent(display_,surface_,surface_,context_))
                throw std::runtime_error("cannot create the EGL OpenGL 3.3 context");
        } catch (...) {release();throw;}
#else
        (void)w;(void)h;throw std::runtime_error("this build has no EGL offscreen support");
#endif
    }
    ~Offscreen() {
#if SPP_MOON_EGL
        release();
#endif
    }
    Offscreen(const Offscreen&)=delete;Offscreen& operator=(const Offscreen&)=delete;
};
}  // namespace spp::app
