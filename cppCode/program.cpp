////////////////////////////////////
#include "program.h"

#include "star.h"
#include "centroid.h"
#include "calccatalog.h"
#include "calckinematic.h"
#include "configloader.h"

#include <list>
#include <vector>
#include <string>

#include <TException.h>
#include <TVersion.h>
#include <TLogger.h>
using namespace TRiOLD;

////////////////////////////////////
void Program::_setup()
{
    m_prgname = "Galaxy Kinematic";
    m_prgversion = Version(1, 2, 0, Version::BETA).toString(true);
    m_argvs = Argvs();
}

Program::~Program()
{

}

Program::Program(int argc, char** argv)
{
    _setup();
    _parseArguments(argc, argv);
}

void Program::_parseArguments(int argc, char** argv)
{
    std::vector<std::string> arguments(argv, argv + argc);
    if(arguments.size() == 1)
    {
        m_argvs.processType = ProcessType::HELP;
        return;
    }

    for(unsigned int i = 0; i < arguments.size(); i++)
    {
        if(arguments.at(i).compare("--version") == 0)
        {
            m_argvs.processType = ProcessType::VERSION;
            continue;
        }
        if(arguments.at(i).compare("--help") == 0)
        {
            m_argvs.processType = ProcessType::HELP;
            continue;
        }
        if(arguments.at(i).compare("--crtstructcatalog") == 0 ||
           arguments.at(i).compare("--catalog") == 0)
        {
            m_argvs.processType = ProcessType::CRTSTRUCTCATALOG;
            continue;
        }
        if(arguments.at(i).compare("--crtpixelstructcatalog") == 0 ||
           arguments.at(i).compare("--pixcatalog") == 0)
        {
            m_argvs.processType = ProcessType::CRTPXLSTRUCTCATALOG;
            continue;
        }
        if(arguments.at(i).compare("--calckinematic") == 0 ||
           arguments.at(i).compare("--kinematic") == 0)
        {
            m_argvs.processType = ProcessType::CALCKINEMATIC;
            continue;
        }

        if(arguments.at(i).compare("--infile") == 0 ||
           arguments.at(i).compare("-i") == 0)
            if(i+1 < arguments.size())
            {
                m_argvs.inFilePath = arguments.at(++i);
                continue;
            }
        if(arguments.at(i).compare("--onfile") == 0 ||
           arguments.at(i).compare("-o") == 0)
            if(i+1 < arguments.size())
            {
                m_argvs.onFilePath = arguments.at(++i);
                continue;
            }

        if(arguments.at(i).compare("--config") == 0 ||
           arguments.at(i).compare("--conf") == 0 ||
           arguments.at(i).compare("-c") == 0)
            if(i+1 < arguments.size())
            {
                m_argvs.configPath = arguments.at(++i);
                continue;
            }

        if(arguments.at(i).compare("--agrstrct") == 0)
            if(i+1 < arguments.size())
            {
                std::string value = arguments.at(++i);
                if(!value.empty())
                   m_argvs.inFileStructIsAgreed = (value.compare("true") == 0 ||
                                                   std::atoi(value.c_str()) != 0);
                continue;
            }

        if(arguments.at(i).compare("--log") == 0 ||
           arguments.at(i).compare("-l") == 0)
            if(i+1 < arguments.size())
            {
                m_argvs.logFilePath = arguments.at(++i);
                continue;
            }

        if(arguments.at(i).compare("--clog") == 0)
            if(i+1 < arguments.size())
            {
                std::string value = arguments.at(++i);
                if(!value.empty())
                   m_argvs.isWithConsoleLog = (value.compare("true") == 0 ||
                                               std::atoi(value.c_str()) != 0);
                continue;
            }
    }
}

void Program::_initLogger()
{
    LOG.setLevel(Logger::Level::INFO);
    LOG.setFlag_isWithCout(m_argvs.isWithConsoleLog);
    if(!m_argvs.logFilePath.empty())
        LOG.setFilepath(m_argvs.logFilePath);
}

void Program::_process_version()
{
    std::cout << m_prgname << std::endl;
    std::cout << m_prgversion << std::endl;
}

void Program::_process_help()
{
    std::cout << m_prgname << std::endl;
    std::cout << m_prgversion << std::endl;
    std::cout << std::endl;

    std::cout <<
        "Process type:" << std::endl <<
        " --version        - show version" << std::endl <<
        " --help           - show help" << std::endl <<
        " --catalog        - create struct catalog" << std::endl <<
        " --pixcatalog     - create struct pixels catalog" << std::endl <<
        " --kinematic      - calculate kinematic parameters" << std::endl <<
        ""  << std::endl <<
        "Files path:" << std::endl <<
        " --infile or -i   - path to input file" << std::endl <<
        " --onfile or -o   - path to onput file" << std::endl <<
        " --config or -c   - path to config file" << std::endl <<
        " --log or -l      - path to log file" << std::endl <<
        ""  << std::endl <<
        "Other variables:" << std::endl <<
        " --agrstrct       - use agreed struct infile? \"true\" or \"false\"" << std::endl <<
        " --clog           - use cout log? \"true\" or \"false\"" << std::endl <<
        ""  << std::endl <<
        "Example:" << std::endl <<
        "program --kinematic -i /path/to/catalog.csv "
        "-o /path/to/kinematic.csv --agrstrct true -clog true" << std::endl << std::endl;
}

void Program::_process_crtStructCatalog()
{
    LOG.writeInfo("Start process: create struct catalog.");

    LOG.writeInfo("Start reading ConfigTable...");
    LOG.writeInfo("FilePath: " + m_argvs.configPath);
    CalcCatalog::ConfigTable config;
    ConfigLoader::load(config, m_argvs.configPath);
    LOG.writeInfo("ConfigTable has been loaded.");

    LOG.writeInfo("Start reading input star catalog...");
    LOG.writeInfo("FilePath: " + m_argvs.inFilePath);
    std::list<Star> stars = CalcCatalog::readCatalog(m_argvs.inFilePath, config);
    LOG.writeInfo("Star catalog has been loaded. Stars amount = " +
                  std::to_string(stars.size()));

    LOG.writeInfo("Start writing result struct star catalog...");
    LOG.writeInfo("FilePath: " + m_argvs.onFilePath);
    CalcCatalog::writeCatalog_agreedStruct(m_argvs.onFilePath, stars);
    LOG.writeInfo("Result struct star catalog has been written.");

    LOG.writeInfo("The process is completed.");
}

void Program::_process_crtPxlStructCatalog()
{
    LOG.writeInfo("Start process: calculate kinematic parameters.");

    std::list<Star> stars;
    if(m_argvs.inFileStructIsAgreed)
    {
        LOG.writeInfo("Start reading input star catalog (agreed struct)...");
        LOG.writeInfo("FilePath: " + m_argvs.inFilePath);
        stars = CalcCatalog::readCatalog_agreedStruct(m_argvs.inFilePath);
    }
    else
    {
        LOG.writeInfo("Start reading ConfigTable...");
        LOG.writeInfo("FilePath: " + m_argvs.configPath);
        CalcCatalog::ConfigTable configTable;
        ConfigLoader::load(configTable, m_argvs.configPath);
        LOG.writeInfo("ConfigTable has been loaded.");

        LOG.writeInfo("Start reading input star catalog...");
        LOG.writeInfo("FilePath: " + m_argvs.inFilePath);
        stars = CalcCatalog::readCatalog(m_argvs.inFilePath, configTable);
    }
    LOG.writeInfo("Star catalog has been loaded. "
                  "Stars amount = " +std::to_string(stars.size()));

    LOG.writeInfo("Start reading ConfigPixelization...");
    LOG.writeInfo("FilePath: " + m_argvs.configPath);
    CalcCatalog::ConfigPixelization configPix;
    ConfigLoader::load(configPix, m_argvs.configPath);
    LOG.writeInfo("ConfigPixelization has been loaded.");

    LOG.writeInfo("Start star catalog pixelization...");
    std::list<Star> pixels = CalcCatalog::createPixCatalog_agreedStruct(stars, configPix);
    std::string pxLen = std::to_string(pixels.size());
    LOG.writeInfo("Star catalog has been pixelated. There are " + pxLen + " pixels.");
    LOG.writeInfo("Start writing result pixelated struct star catalog...");
    LOG.writeInfo("FilePath: " + m_argvs.onFilePath);
    CalcCatalog::writeCatalog_agreedStruct(m_argvs.onFilePath, pixels);
    LOG.writeInfo("Result pixelated struct star catalog has been written.");

    LOG.writeInfo("The process is completed.");
}

void Program::_process_calcKinematic()
{
    LOG.writeInfo("Start process: calculate kinematic parameters.");

    std::list<Star> stars;
    if(m_argvs.inFileStructIsAgreed)
    {
        LOG.writeInfo("Start reading input star catalog (agreed struct)...");
        LOG.writeInfo("FilePath: " + m_argvs.inFilePath);
        stars = CalcCatalog::readCatalog_agreedStruct(m_argvs.inFilePath);
    }
    else
    {
        LOG.writeInfo("Start reading ConfigTable...");
        LOG.writeInfo("FilePath: " + m_argvs.configPath);
        CalcCatalog::ConfigTable configTable;
        ConfigLoader::load(configTable, m_argvs.configPath);
        LOG.writeInfo("ConfigTable has been loaded.");

        LOG.writeInfo("Start reading input star catalog...");
        LOG.writeInfo("FilePath: " + m_argvs.inFilePath);
        stars = CalcCatalog::readCatalog(m_argvs.inFilePath, configTable);
    }
    LOG.writeInfo("Star catalog has been loaded. "
                  "Stars amount = " +std::to_string(stars.size()));


    LOG.writeInfo("Start reading ConfigConstants...");
    LOG.writeInfo("FilePath: " + m_argvs.configPath);
    CalcCatalog::ConfigConstants configConst;
    ConfigLoader::load(configConst, m_argvs.configPath);
    LOG.writeInfo("ConfigConstants has been loaded.");

    LOG.writeInfo("Start reading ConfigProcessing...");
    LOG.writeInfo("FilePath: " + m_argvs.configPath);
    CalcKinematic::ConfigProcessing configProc;
    ConfigLoader::load(configProc, m_argvs.configPath);
    LOG.writeInfo("ConfigProcessing has been loaded.");

    LOG.writeInfo("Start calculating kinematic parameters for centroids...");
    std::vector<Centroid> centriods = CalcKinematic::calcCentroids(stars, configProc);
    LOG.writeInfo("Kinematic parameters has been calculated.. "
                  "Centriods amount = " +std::to_string(centriods.size()));

    LOG.writeInfo("Start writing result kinematic parameters...");
    LOG.writeInfo("FilePath: " + m_argvs.onFilePath);
    CalcCatalog::writeCentroids(m_argvs.onFilePath, centriods, configConst);
    LOG.writeInfo("Result kinematic parameters has been written.");

    LOG.writeInfo("The process is completed.");
}

int Program::process()
{
    try
    {
        switch (m_argvs.processType)
        {
            case Program::ProcessType::VERSION:
                _process_version();
            break;

            case Program::ProcessType::HELP:
                _process_help();
            break;

            case Program::ProcessType::CRTSTRUCTCATALOG:
                _initLogger();
                _process_crtStructCatalog();
            break;

            case Program::ProcessType::CRTPXLSTRUCTCATALOG:
                _initLogger();
                _process_crtPxlStructCatalog();
            break;

            case Program::ProcessType::CALCKINEMATIC:
                _initLogger();
                _process_calcKinematic();
            break;

            default: throw Exception("Unknown process type!", -10);
        }
    }
    catch (const Exception& exc)
    {
        LOG.writeFatal(exc.getComment());
        return exc.getCode();
    }
    catch (const std::exception& exc)
    {
        LOG.writeFatal(exc.what());
        return -1;
    }
    catch (...)
    {
        LOG.writeFatal("UNKNOWN ERROR");
        return -1;
    }
    return 0;
}

////////////////////////////////////
