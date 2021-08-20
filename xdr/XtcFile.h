#include "XdrWrapper.h"
#include "libxdr/xtc_seek.h"

#include <chrono>

class XtcFile:public XdrWrapper
{
    public:
        XtcFile(const XdrInput& input);
        virtual ~XtcFile(){};

        virtual bool readFrame(int FrameNum) override;
        virtual void readNumAtoms() override;
        virtual void readNframes() override;
    private:
        Frame::Real precision_;

        std::vector<int64_t> offsets_;
};