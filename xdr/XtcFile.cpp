#include "XtcFile.h"
namespace XdrFiles
{
    registry_<XtcFile> register_Xtc("xtc");
}

XtcFile::XtcFile(const XdrInput& input)
:XdrWrapper(input)
{}

void XtcFile::readNumAtoms()
{
    ASSERT((isOpen()), "The file is not opened.");

    int success = read_xtc_natoms(const_cast<char*>(path_.c_str()), &natoms_);

    ASSERT((success == exdrOK), "The process to read xtc natoms is not sucessful.");
}

void XtcFile::readNframes()
{
    int est_nframes=0;
    int64_t* offsets=nullptr;

    read_xtc_n_frames(const_cast<char*>(path_.c_str()),&nframes_, &est_nframes, &offsets);
    offsets_.insert(offsets_.end(),offsets, offsets+nframes_);
}

bool XtcFile::readFrame(int FrameNum)
{
    ASSERT((isOpen()), "The file is not opened.");
    if (FrameNum > nframes_)
    {
        return true;
    }

    auto& positions_ = frame_.accessPositions();
    ASSERT((positions_.size() == natoms_), "The shape of positions is not of natoms size.");

    int step;
    Frame::Real time;
    matrix box_;
    Frame::Matrix matrix_box_;

    auto positions_ptr_ = (rvec*)positions_.data();

    xdr_seek(file_,offsets_[FrameNum], 0);
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