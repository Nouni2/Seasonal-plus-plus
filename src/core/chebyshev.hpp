#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include "spp/core/constants.hpp"

namespace spp::detail {

// Offset Chebyshev interpolation of a fixed-width vector on [-1,1]. Keeping
// the central sample separate avoids fitting an astronomical coordinate's
// large constant part into every coefficient. The caller owns the domain.
template<std::size_t Width,std::size_t Degree>
struct ChebyshevSeries {
    static constexpr std::size_t nodes=Degree+1;
    using Value=std::array<double,Width>;
    Value offset{};
    std::array<Value,nodes> coefficients{};

    template<class Sample>
    static ChebyshevSeries fit(Sample sample) {
        struct Projection {
            std::array<double,nodes> x;
            std::array<std::array<double,nodes>,nodes> weight;
            Projection() {
                for(std::size_t k=0;k<nodes;++k) {
                    x[k]=std::cos(pi*(k+0.5)/nodes);
                    for(std::size_t j=0;j<nodes;++j)weight[j][k]=(2.0/nodes)*std::cos(pi*j*(k+0.5)/nodes);
                }
            }
        };
        static const Projection projection;
        std::array<Value,nodes> values;
        for(std::size_t k=0;k<nodes;++k)values[k]=sample(projection.x[k]);
        ChebyshevSeries result;result.offset=values[nodes/2];
        for(std::size_t j=0;j<nodes;++j)for(std::size_t k=0;k<nodes;++k)
            for(std::size_t q=0;q<Width;++q)
                result.coefficients[j][q]+=projection.weight[j][k]*(values[k][q]-result.offset[q]);
        for(auto& c:result.coefficients[0])c*=0.5;
        return result;
    }
    Value operator()(double x)const {
        Value b1{},b2{};
        for(std::size_t j=nodes;j-->1;)for(std::size_t q=0;q<Width;++q) {
            const double b=coefficients[j][q]+2*x*b1[q]-b2[q];b2[q]=b1[q];b1[q]=b;
        }
        Value result;
        for(std::size_t q=0;q<Width;++q)result[q]=offset[q]+(coefficients[0][q]+x*b1[q]-b2[q]);
        return result;
    }
};
} // namespace spp::detail
