#include <assert.h>
#include <fstream>
#include <iostream>
#include <vector>

#include "MeshSim.h"
#include "MeshTypes.h"
#include "ModelEnums.h"
#include "ModelTypes.h"
#include "SimAdvModel.h"
#include "SimInfo.h"
#include "SimModel.h"
#include "SimPList.h"
#include "SimUtil.h"


void messageHandler(int type, const char* msg);


int main(int argc, char* argv[])
{
    if (argc != 4)
    {
        std::cerr << "Usage: " << argv[0]
            << " <model_name> <mesh_name> <wall file name>\n";
        return 1;
    }


    const char* modelName = argv[1];
    const char* meshName = argv[2];
    const char* wallFile = argv[3];


    try
    {
        Sim_logOn("simmeshWboundaryNodes.log");
        SimModel_start(); // Call before Sim_readLicenseFile
        // NOTE: Sim_readLicenseFile() is for internal testing only.  To use,
        // pass in the location of a file containing your keys.  For a release
        // product, use Sim_registerKey()
        Sim_readLicenseFile(0);
        // Tessellation of GeomSim geometry requires Meshing to have started
        MS_init();
        Sim_setMessageHandler(messageHandler);
        pProgress progress = Progress_new();
        Progress_setDefaultCallback(progress);

        pGModel orig_model = GM_load(modelName, NULL, progress);
        std::cout << "Model loaded: " << modelName << std::endl;










        // ************* Clean Up *****************//
        GM_release(orig_model);
        Progress_delete(progress);
        MS_exit();
        Sim_unregisterAllKeys();
        SimModel_stop();
        Sim_logOff();
    }
    catch (pSimInfo err)
    {
        std::cerr << "SimModSuite error caught:" << std::endl;
        std::cerr << "  Error code: " << SimInfo_code(err) << std::endl;
        std::cerr << "  Error string: " << SimInfo_toString(err) << std::endl;
        SimInfo_delete(err);
        return 1;
    }
    catch
    (...)
    {
        std::cerr << "Unhandled exception caught" << std::endl;
        return 1;
    }


    return 0;
}


void messageHandler(int type, const char* msg)
{
    switch (type)
    {
    case Sim_InfoMsg:
        std::cout << "Info: " << msg << std::endl;
        break;
    case Sim_DebugMsg:
        std::cout << "Debug: " << msg << std::endl;
        break;
    case Sim_WarningMsg:
        std::cout << "Warning: " << msg << std::endl;
        break;
    case Sim_ErrorMsg:
        std::cout << "Error: " << msg << std::endl;
        break;
    }
}
