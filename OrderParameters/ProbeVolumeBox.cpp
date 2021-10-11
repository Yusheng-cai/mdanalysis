#include "ProbeVolumeBox.h"

namespace ProbeVolumes
{
    static const registry_<ProbeVolumeBox> register_box("box");
}

ProbeVolumeBox::ProbeVolumeBox(ProbeVolumeInput& input)
:ProbeVolume(input)
{
    input.ParamPack.ReadArrayNumber("xrange", ParameterPack::KeyType::Required, xrange_);
    input.ParamPack.ReadArrayNumber("yrange", ParameterPack::KeyType::Required, yrange_);
    input.ParamPack.ReadArrayNumber("zrange", ParameterPack::KeyType::Required, zrange_);
    input.ParamPack.ReadNumber("sigma", ParameterPack::KeyType::Optional, sigma_);
    input.ParamPack.ReadNumber("alphac", ParameterPack::KeyType::Optional, ac_);


    ASSERT((! isDynamic()), "The box probe volume does not yet support dynamic atom group.");

    setGeometry(); 
}

void ProbeVolumeBox::setGeometry()
{   
    // Initialize the indicator functions
    dx_ = xrange_[1] - xrange_[0];
    dy_ = yrange_[1] - yrange_[0];
    dz_ = zrange_[1] - zrange_[0];

    center_[0] = 0.5*(xrange_[0] + xrange_[1]);
    center_[1] = 0.5*(yrange_[0] + yrange_[1]);
    center_[2] = 0.5*(zrange_[0] + zrange_[1]);

    // Set indicator function such that the center of the PV is at (0,0,0)
    func_[0] = IndicatorFunction2d(sigma_, ac_, -dx_/2, dx_/2);
    func_[1] = IndicatorFunction2d(sigma_, ac_, -dy_/2, dy_/2);
    func_[2] = IndicatorFunction2d(sigma_, ac_, -dz_/2, dz_/2);
}

ProbeVolumeOutput ProbeVolumeBox::calculate(const Real3& x) const
{
    // First calculate the pbc corrected distance of x with respect to the center of the Box PV
    Real3 dist;
    Real sq_dist;
    simbox_.calculateDistance(x, center_, dist, sq_dist); 

    // Initialize the output variable
    ProbeVolumeOutput output;
    Real3 htilde;
    Real3 h;
    Real3 dhtilde;

    // Find h(x), h(y) , h(z) && h'(x) , h'(y) , h'(z)
    for (int i=0;i<3;i++)
    {
        Real htilde_x, h_x, dhtilde_dx;
        func_[i].calculate(dist[i],h_x, htilde_x, dhtilde_dx);        
        htilde[i] = htilde_x;
        dhtilde[i] = dhtilde_dx;
        h[i] = h_x;
    }

    // Calculate the derivative dh/dx = h'(x) * h(y) * h(z)
    dhtilde[0] = dhtilde[0] * htilde[1] * htilde[2];
    dhtilde[1] = dhtilde[1] * htilde[0] * htilde[2];
    dhtilde[2] = dhtilde[2] * htilde[0] * htilde[1];

    output.dhtilde_dx_ = dhtilde;
    output.htilde_x_ = htilde[0] * htilde[1] * htilde[2];
    output.hx_ = h[0] * h[1] * h[2];

    return output;
}