#include "XdrWrapper.h"
XdrWrapper::XdrWrapper(const XdrInput& input)
:pack_(input.pack), apath_(input.apath_)
{
    pack_.ReadString("path", ParameterPack::KeyType::Required, path_);
    path_ = FileSystem::joinPath(apath_, path_);

    std::cout << "path = " << path_ << std::endl;

    bool readmode = pack_.ReadString("mode", ParameterPack::KeyType::Optional, operation_mode_);
 

    if ( ! readmode)
    {
        operation_mode_ = "read";
    }
}

void XdrWrapper::open()
{
    std::string mode_;

    ASSERT((operation_mode_ == "read" || operation_mode_ == "append" || operation_mode_ == "write"), "mode has to be one of 'read' & 'append' & 'write'");

    // default is reading mode
    if (operation_mode_ == "read")
    {
        mode_ = "r";
    }
    else if (operation_mode_ == "append")
    {
        mode_ = "a";
    }   
    else if (operation_mode_ == "write")
    {
        mode_ = "w";
    }

    file_ = xdrfile_open(path_.c_str(), mode_.c_str());
    ASSERT((isOpen()), "The file did not open correctly.");

    // Read Number of Atoms 
    readNumAtoms();

    frame_.setNumAtoms(natoms_);

    int nframes = 0;
    int est_nframes = 0;
    int64_t* offsets = nullptr;

    // readNframes();
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