#include "XdrWrapper.h"

void XdrWrapper::open(std::string fname, Mode mode)
{
    name_ = fname;

    switch(mode)
    {
        case Mode::Read:
            operation_mode_="r";
            break;
        case Mode::Write:
            operation_mode_="w";
            break;
        case Mode::Append:
            operation_mode_="a";
            break;
    };

    file_ = xdrfile_open(fname.c_str(), operation_mode_.c_str());
    ASSERT((isOpen()), "The file did not open correctly.");

    // Read Number of Atoms 
    readNumAtoms();

    frame_.setNumAtoms(natoms_);
}

bool XdrWrapper::isOpen()
{
    if (file_ != nullptr)
    {
        return true;
    }
    else
    {
        return false;
    }
}

void XdrWrapper::close()
{
    ASSERT((isOpen()), "You are trying to close the file while it wasn't opened in the first place.");
    int ret = xdrfile_close(file_);

    ASSERT((ret == 0), "The xdr file was closed incorrectly.");
}

XdrWrapper::~XdrWrapper()
{
    close();
}