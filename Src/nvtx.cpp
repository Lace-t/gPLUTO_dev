// #ifdef _USE_NVTX// SPDX-License-Identifier: BSD-3-Clause
#ifdef _USE_NVTX
#include <nvtx3/nvToolsExt.h>
#include <string.h>
#endif

#ifdef _USE_ROCTX
//#include <rocprofiler-sdk-roctx/roctx.h>
#include <roctracer/roctx.h>
#include <string.h>
#endif

void mynvtxstart_(const char *name, int color) {
#ifdef _USE_NVTX
   nvtxEventAttributes_t eventAttrib = {0};
   eventAttrib.version = NVTX_VERSION;
   eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE;
   eventAttrib.colorType = NVTX_COLOR_ARGB;
   eventAttrib.color = color;
   eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII;
   eventAttrib.message.ascii = name;
   nvtxRangePushEx(&eventAttrib);
#elif defined(_USE_ROCTX)
    //printf("ROCTX: %s\n", name);
    /*roctxEventAttributes_t eventAttrib = {0};
    eventAttrib.version = ROCTX_VERSION;
    eventAttrib.size = ROCTX_EVENT_ATTRIB_STRUCT_SIZE;
    eventAttrib.colorType = ROCTX_COLOR_ARGB;
    eventAttrib.color = color;
    eventAttrib.messageType = ROCTX_MESSAGE_TYPE_ASCII;
    eventAttrib.message.ascii = name;*/
    roctxRangePush(name);
#endif
}

void mynvtxstop_() {
#ifdef _USE_NVTX
    nvtxRangePop();
#elif defined(_USE_ROCTX)
    roctxRangePop();
#endif
}
