////////////////////////////////////
#ifndef CONFIGLOADER_H
#define CONFIGLOADER_H

////////////////////////////////////
#include <string>

#include "calccatalog.h"
#include "calckinematic.h"

////////////////////////////////////
namespace ConfigLoader
{
    bool load(CalcCatalog::ConfigTable &config, const std::string &xmlpath);
    bool load(CalcCatalog::ConfigConstants &config, const std::string &xmlpath);
    bool load(CalcCatalog::ConfigPixelization &config, const std::string &xmlpath);
    bool load(CalcKinematic::ConfigProcessing &config, const std::string &xmlpath);
};

////////////////////////////////////
#endif // CONFIGLOADER_H
