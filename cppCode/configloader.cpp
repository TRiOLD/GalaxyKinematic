////////////////////////////////////
#include "configloader.h"

//pugixml Lib
#include <pugixml.hpp>
using namespace pugi;

////////////////////////////////////
bool _check(const xml_attribute& att, const std::string& value = "")
{
    if(!att)
        return false;
    if(std::string(att.value()).empty())
        return false;
    if(!value.empty())
        if(value.compare(att.value()) != 0)
            return false;
    return true;
}

bool _check(const xml_node& node, const std::string& type, const std::string& dim = "")
{
    if(!node)
        return false;

    if(!_check(node.attribute("value")))
        return false;

    if(!_check(node.attribute("type"), type))
        return false;
    if(!dim.empty())
        if(!_check(node.attribute("dim"), dim))
            return false;
    return true;
}

bool ConfigLoader::load(CalcCatalog::ConfigTable& config, const std::string& xmlpath)
{
    xml_document doc;
    if (!doc.load_file(xmlpath.c_str()))
        return false;

    xml_node configNode = doc.child("GalaxyKinematic").child("Catalog").child("ConfigTable");
    if(!configNode)
        return false;

    xml_node param;
    param = configNode.child("fileType");
    if(_check(param, "int"))
        config.fileType = param.attribute("value").as_int();

    param = configNode.child("isWithHeader");
    if(_check(param, "bool"))
        config.isWithHeader = param.attribute("value").as_bool();

    param = configNode.child("coordInRad");
    if(_check(param, "bool"))
        config.coordInRad = param.attribute("value").as_bool();

    param = configNode.child("cn_ra");
    if(_check(param, "int"))
        config.cn_ra = param.attribute("value").as_int();

    param = configNode.child("cn_de");
    if(_check(param, "int"))
        config.cn_de = param.attribute("value").as_int();

    param = configNode.child("cn_plx");
    if(_check(param, "int"))
        config.cn_plx = param.attribute("value").as_int();

    param = configNode.child("cn_muRa");
    if(_check(param, "int"))
        config.cn_muRa = param.attribute("value").as_int();

    param = configNode.child("cn_muDe");
    if(_check(param, "int"))
        config.cn_muDe = param.attribute("value").as_int();

    param = configNode.child("cn_Vr");
    if(_check(param, "int"))
        config.cn_Vr = param.attribute("value").as_int();

    return true;
}

bool ConfigLoader::load(CalcCatalog::ConfigPixelization &config, const std::string &xmlpath)
{
    xml_document doc;
    if (!doc.load_file(xmlpath.c_str()))
        return false;

    xml_node configNode = doc.child("GalaxyKinematic").child("Catalog").child("Pixelization");
    if(!configNode)
        return false;

    xml_node param;
    param = configNode.child("minX");
    if(_check(param, "double", "kpc"))
        config.minX = param.attribute("value").as_double();

    param = configNode.child("maxX");
    if(_check(param, "double", "kpc"))
        config.maxX = param.attribute("value").as_double();

    param = configNode.child("minY");
    if(_check(param, "double", "kpc"))
        config.minY = param.attribute("value").as_double();

    param = configNode.child("maxY");
    if(_check(param, "double", "kpc"))
        config.maxY = param.attribute("value").as_double();

    param = configNode.child("minZ");
    if(_check(param, "double", "kpc"))
        config.minZ = param.attribute("value").as_double();

    param = configNode.child("maxZ");
    if(_check(param, "double", "kpc"))
        config.maxZ = param.attribute("value").as_double();

    param = configNode.child("pixelHalfWidth");
    if(_check(param, "double", "kpc"))
        config.pixelHalfWidth = param.attribute("value").as_double();

    param = configNode.child("minStarCount");
    if(_check(param, "int"))
        config.minStarCount = param.attribute("value").as_int();

    return true;
}

bool ConfigLoader::load(CalcCatalog::ConfigConstants& config, const std::string& xmlpath)
{
    xml_document doc;
    if (!doc.load_file(xmlpath.c_str()))
        return false;

    xml_node configNode = doc.child("GalaxyKinematic").child("Catalog").child("ConfigConstants");
    if(!configNode)
        return false;

    xml_node param;
    param = configNode.child("R_Sun");
    if(_check(param, "double", "kpc"))
        config.R_Sun = param.attribute("value").as_double();

    param = configNode.child("VX_Sun");
    if(_check(param, "double", "km/sec"))
        config.VX_Sun = param.attribute("value").as_double();

    param = configNode.child("VY_Sun");
    if(_check(param, "double", "km/sec"))
        config.VY_Sun = param.attribute("value").as_double();

    param = configNode.child("VZ_Sun");
    if(_check(param, "double", "km/sec"))
        config.VZ_Sun = param.attribute("value").as_double();

    return true;
}

bool ConfigLoader::load(CalcKinematic::ConfigProcessing& config, const std::string& xmlpath)
{
    xml_document doc;
    if (!doc.load_file(xmlpath.c_str()))
        return false;

    xml_node configNode = doc.child("GalaxyKinematic").child("Kinematic").child("ConfigProcessing");
    if(!configNode)
        return false;

    xml_node param;
    param = configNode.child("threadsAmount");
    if(_check(param, "int"))
        config.threadsAmount = param.attribute("value").as_int();

    param = configNode.child("RAMlimit");
    if(_check(param, "int", "MB"))
        config.RAMlimit = param.attribute("value").as_int();
    if(_check(param, "int", "GB"))
        config.RAMlimit = param.attribute("value").as_int()*1024;

    param = configNode.child("minX");
    if(_check(param, "double", "kpc"))
        config.minX = param.attribute("value").as_double();

    param = configNode.child("maxX");
    if(_check(param, "double", "kpc"))
        config.maxX = param.attribute("value").as_double();

    param = configNode.child("minY");
    if(_check(param, "double", "kpc"))
        config.minY = param.attribute("value").as_double();

    param = configNode.child("maxY");
    if(_check(param, "double", "kpc"))
        config.maxY = param.attribute("value").as_double();

    param = configNode.child("minZ");
    if(_check(param, "double", "kpc"))
        config.minZ = param.attribute("value").as_double();

    param = configNode.child("maxZ");
    if(_check(param, "double", "kpc"))
        config.maxZ = param.attribute("value").as_double();

    param = configNode.child("step");
    if(_check(param, "double", "kpc"))
        config.step = param.attribute("value").as_double();

    param = configNode.child("starsRegionRadius");
    if(_check(param, "double", "kpc"))
        config.starsRegionRadius = param.attribute("value").as_double();

    return true;
}

////////////////////////////////////
