#include "XtcFile.h"
namespace XdrFiles
{
    registry_<XtcFile> register_Xtc("xtc");
}

XtcFile::XtcFile()
:XdrWrapper()
{}

void XtcFile::readNumAtoms()
{
    ASSERT((isOpen()), "The file is not opened.");

    int success = read_xtc_natoms(const_cast<char*>(name_.c_str()), &natoms_);

    ASSERT((success == exdrOK), "The process to read xtc natoms is not sucessful.");
}

void XtcFile::readNframes()
{
    int est_nframes=0;
    int64_t* offsets=nullptr;

    read_xtc_n_frames(const_cast<char*>(name_.c_str()),&nframes_, &est_nframes, &offsets);
}

bool XtcFile::readNextFrame()
{
    ASSERT((isOpen()), "The file is not opened.");

    auto& positions_ = frame_.accessPositions();
    ASSERT((positions_.size() == natoms_), "The shape of positions is not of natoms size.");

    int step;
    Frame::Real time;
    matrix box_;
    Frame::Matrix matrix_box_;

    auto positions_ptr_ = (rvec*)positions_.data();
    int ret = read_xtc(file_, natoms_, &step, &time, box_, positions_ptr_, &precision_);
    ASSERT((ret == exdrOK || ret == exdrENDOFFILE), "The reading operation in xtc file is not sucessful.");

    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            matrix_box_[i][j] = box_[i][j];
        }
    }
    frame_.setBox(matrix_box_);
    frame_.setTime(time);
    frame_.setStep(step);

    if (ret == exdrENDOFFILE)
    {
        return true;
    }
    else
    {
        return false;
    }
}