#include "BetaFactorWriter.h"

BetaFactorWriter::BetaFactorWriter(ParameterPack& input)
{
    input.ReadString("name", ParameterPack::KeyType::Required, filename_);

    ofs_.open(filename_);
}

void BetaFactorWriter::write(int frameNum, const std::vector<Real>& data)
{
    ofs_ << frameNum << " ";

    for (int i=0;i<data.size();i++)
    {
        ofs_ << data[i] << " ";
    }

    ofs_ << "\n";
}