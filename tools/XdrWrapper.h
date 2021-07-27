#pragma once
#include "libxdr/xdrfile_trr.h"
#include "libxdr/xdrfile_xtc.h"
#include "libxdr/xdrfile.h"
#include "Assert.h"
#include "Frame.h"
#include "GenericFactory.h"

#include <string>

class XdrWrapper
{
    public:
        enum Mode
        {
            Read, Write, Append
        };

        XdrWrapper(){};    
        void open(std::string name, Mode);
        virtual ~XdrWrapper();

        // close the file
        void close();

        // Check if the xdr file is opened
        bool isOpen();

        int getNumAtoms() const {return natoms_;}

        virtual void readNextFrame() = 0;
        virtual void readNumAtoms()  = 0;

    protected:
        XDRFILE* file_=nullptr;
        int natoms_;
        std::string name_;
        std::string operation_mode_;

        Frame frame_;
};

namespace XdrFiles
{
    using key  = std::string;
    using Base = XdrWrapper;

    using factory = GenericFactory<Base, key>;

    template<class D>
    using registry_ = RegisterInFactory<Base, D, key>;
};