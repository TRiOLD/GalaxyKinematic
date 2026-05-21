# GalaxyKinematic

## Overview

GalaxyKinematic is a program for processing astronomical data, calculating kinematic parameters, and generating structured catalogs. The project relies on two external repositories:

- **CppHelperClasses** - [GitHub Repository](https://github.com/TRiOLD/CppHelperClasses)
- **CppExternalSourceLibs** - [GitHub Repository](https://github.com/TRiOLD/CppExternalSourceLibs)

## Project Structure

```
CppExternalSourceLibs/
├── alglib/
│   ├── 3.19.0/
│   │   ├── include/ ...
│   │   ├── source/ ...
│   │   ├── CMakeLists.txt
├── pugixml/
│   ├── 1.14/
│   │   ├── include/ ...
│   │   ├── source/ ...
│   │   ├── CMakeLists.txt
CppHelperClasses/ ...
GalaxyKinematic/
│── cppCode/ ...
│── xml/ ...
│── CMakeLists.txt
```

## Installation

1. Clone the repositories:
   ```sh
   git clone https://github.com/TRiOLD/GalaxyKinematic.git
   git clone https://github.com/TRiOLD/CppHelperClasses.git
   git clone https://github.com/TRiOLD/CppExternalSourceLibs.git
   ```
2. Build the project using CMake:
   ```sh
   mkdir build && cd build
   cmake ..
   cmake --build .
   ```

## Usage

GalaxyKinematic supports several processing modes and configurable file paths.

### **Process Type**

```
--version        - Show version
--help           - Show help
--catalog        - Create structured catalog
--pixcatalog     - Create structured pixels catalog
--kinematic      - Calculate kinematic parameters
```

### **File Paths**

```
--infile or -i   - Path to input file
--onfile or -o   - Path to output file
--config or -c   - Path to config file
--log or -l      - Path to log file
```

### **Other Variables**

```
--agrstrct       - Use agreed struct infile? "true" or "false"
--clog           - Use cout log? "true" or "false"
```

### **Example Usage**

Linux:
```sh
./GalaxyKinematic --kinematic -i /path/to/catalog.csv -o /path/to/kinematic.csv --agrstrct true --clog true
```
Windows:
```sh
GalaxyKinematic.exe --kinematic -i /path/to/catalog.csv -o /path/to/kinematic.csv --agrstrct true --clog true
```

## License

Free to use.

