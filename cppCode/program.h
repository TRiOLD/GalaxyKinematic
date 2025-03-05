////////////////////////////////////
#ifndef PROGRAM_H
#define PROGRAM_H

////////////////////////////////////
#include <vector>
#include <string>

////////////////////////////////////
class Program
{
    ////////////////////////////////
public:
    enum ProcessType
    {
        UNKNOWN,
        VERSION,            // --version
        HELP,               // --help
        CRTSTRUCTCATALOG,   // --create struct catalog
        CRTPXLSTRUCTCATALOG,// --create pixel struct catalog
        CALCKINEMATIC       // --calculate kinematic parameters
    };

    struct Argvs
    {
        ProcessType processType = UNKNOWN;
        std::string configPath = "";
        std::string inFilePath = "";
        std::string onFilePath = "";

        std::string logFilePath = "";
        bool isWithConsoleLog = true;

        bool inFileStructIsAgreed = true;
    };

private:
    void _setup();

public:
    ~Program();
    Program(int argc, char** argv);

private:
    std::string m_prgname;
    std::string m_prgversion;
    Argvs m_argvs;

private:
    void _parseArguments(int argc, char** argv);
    void _initLogger();

    void _process_version();
    void _process_help();
    void _process_crtStructCatalog();
    void _process_crtPxlStructCatalog();
    void _process_calcKinematic();

public:
    int process();
};

////////////////////////////////////
#endif // PROGRAM_H
