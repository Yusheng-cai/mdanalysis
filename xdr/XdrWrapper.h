#pragma once
#ifndef XDR_WRAPPER_
#define XDR_WRAPPER_
#include "libxdr/xdrfile_trr.h"
#include "libxdr/xdrfile_xtc.h"
#include "libxdr/xdrfile.h"
#include "tools/Assert.h"
#include "Frame.h"
#include "tools/GenericFactory.h"
#include "tools/CommonTypes.h"

#include <string>
#include <vector>
#include <iostream>


class XdrWrapper
{
    public:
        using Real = CommonTypes::Real;
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

        virtual bool readNextFrame() = 0;
        virtual void readNumAtoms()  = 0;
        const Frame::VectorReal3& getPositions() const{return frame_.getPositions();}
        const Frame::VectorReal3& getVelocities() const{return frame_.getVelocities();}
        const Frame::VectorReal3& getFroces() const{return frame_.getForces();}
        Real getTime() const {return frame_.getTime();}
        int getStep() const {return frame_.getStep();}
        

    protected:
        XDRFILE* file_=nullptr;
        int natoms_;
        std::string name_;
        std::string operation_mode_;

        Frame frame_;
};

namespace XdrFiles
{
    using Key  = std::string;
    using Base = XdrWrapper;

    using factory = GenericFactory<Base, Key>;

    template<class D>
    using registry_ = RegisterInFactory<Base, D, Key>;
};
#endif