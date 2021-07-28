#include "XdrWrapper.h"

class XtcFile:public XdrWrapper
{
    public:
        XtcFile();
        virtual ~XtcFile(){};

        virtual bool readNextFrame() override;
        virtual void readNumAtoms() override;
    private:
        Frame::Real precision_;
};