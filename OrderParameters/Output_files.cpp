#include "Output_files.h"

OutputStream::OutputStream(const ParameterPack& input, const Driver& driver)
:driver_(driver)
{
    ASSERT((input.get_packname() == "output_file"), "The passed in parameter pack is not for Output_file!");

    // Read the name of the file, default to "KMC.out"
    input.ReadString("name", ParameterPack::KeyType::Optional, name_);

    // Read the name of the stride, default to 1
    input.ReadNumber("stride", ParameterPack::KeyType::Optional, stride_);
    ASSERT((stride_ > 0), "Stride must be larger than 0, but provided stride is " << stride_);

    // Read the precision
    input.ReadNumber("precision", ParameterPack::KeyType::Optional, precision_);
    ASSERT((precision_ > 0), "Precision must be larger than 0, but provided stride is " << stride_);


    // Read the value to be outputted
    input.ReadVectorString("value", ParameterPack::KeyType::Required,value_);

    // Open the file with the given name 
    Open();
}

void OutputStream::Open()
{
    // if we are not restarting, then open the file, else we append
    if(! is_restart_)
    {
        ofs_.open(name_);
    }
    else
    {
        ofs_.open(name_, std::ios_base::app);
    }

    ASSERT((isOpen()), "The file with file name " << name_ << " is not found.");

    write_header();
}

void OutputStream::write_header()
{
    ASSERT((isOpen()), "The file with filename " << name_ << " is not found.");

    ofs_ << "#time";

    ofs_ << " ";
    
    for (auto val :value_)
    {
        ofs_ << val;
        ofs_ << "  ";
    }

    ofs_ << "\n";
}

void OutputStream::close()
{
    ASSERT((isOpen()), "The file with filename " << name_ << " is not found.");

    // close the file
    ofs_.close();
}

bool OutputStream::onStep()
{
    int step = driver_.getStep();

    if (step % stride_ == 0)
    {
        return true;
    }
    else
    {
        return false;
    }
}

void OutputStream::printIfOnStep()
{
    // Every simulation must have a time component
    ofs_ << std::fixed << std::setprecision(precision_);

    if (onStep())
    {
        ofs_ << driver_.getTime();

        ofs_ << " ";

        for (auto val:value_)
        {
            // ofs_ << driver_.get_output_val(val);
            
            ofs_ << " ";
        }

        ofs_ << "\n";
    }
}