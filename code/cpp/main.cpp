#include "Viewer.h"

bool parseCmdLine(int argc, char** argv, 
                  ot::Options& opt, 
                  char* & meshFile, 
                  double& scale,
                  int   & verbose)
{
    int arg;
    while ((arg = getopt(argc, argv, "t:i:g:u:a:e:m:v:cf:d:s:")) != -1)
    {
        switch (arg)
        {
            // Options
            case 't':
                opt.tolerance = atof(optarg);
                break;
            case 'i':
                opt.maxIters = atoi(optarg);
                break;
            case 'g':
                opt.gamma = atof(optarg);
                break;
            case 'u':
                opt.upperEntropy = atof(optarg);
                break;
            case 'e':
                opt.epsilon = atof(optarg);
                break;
            case 'd':
                opt.diffIters = atoi(optarg);
                break;
            // Mesh file
            case 'm':
                meshFile = optarg;
                break;
            // Extra parameters
            case 's':
                scale = atof(optarg);
                break;
            case 'v':
                verbose = atoi(optarg);
                break;
        }
    }

    if (meshFile == 0) 
        std::cout << "Usage: " << argv[0] 
        << " -m mesh.obj"
        << " [-t tolerance=1e-6]"
        << " [-i maxIters=1000]"
        << " [-d diffusionIter=10]"
        << " [-g gamma=0]"
        << " [-e epsilon=1e-20]"
        << " [-u upperEntropy=1]"
        << " [-s scale=1]"
        << " [-v verbose=0]"
        << std::endl;

    return meshFile;
}

int main(int argc, char** argv)
{
    ot::Options opt;
    int    verbose  = 0;
    char*  meshFile = 0;
    double scale    = 1.0;
    if (!parseCmdLine(argc, argv, opt, meshFile, scale, verbose)) return 0;

    // load mesh
    TriMesh mesh;
    mesh.read(meshFile);
    mesh.normalize(); // make Area = 1
    std::cout << "MeshArea = " << mesh.computeTotalArea() << std::endl;
    
    // timestep proportinal to mesh size
    double h = mesh.computeMaxEdgeLength();
    std::cout << "h = " << h << std::endl;

    // set gamma a la [Crane et al. 2013]
    if (opt.gamma == 0.) opt.gamma = scale*h*h;
    std::cout << "gamma = " << opt.gamma << std::endl;
 
    // bake kernel
    VectorXd area;
    mesh.computeVertArea(area);

    SparseMatrix matrix;
    mesh.buildHeatKernel(matrix, opt.gamma);

    LinearSolver lsolver;
    lsolver.factorizePosDef(matrix);

    // set OT solver
    ot::ConvSolver otsolver(opt, area, lsolver);

    // gui
    Viewer viewer;
    viewer.meshPtr = &mesh;
    viewer.verbose = verbose;
    viewer.solverPtr = &otsolver;
    viewer.clearData();
    viewer.init(argc,argv);
    return 1;
}
