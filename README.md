## ProjectA

# TODO
- Complete the abPOA.cpp and .hpp files
- Add PaSGAL and vargas

# Compilation
There are two possibilities to compile the code.
To compile the basic code issue the following command in the projectA directory:

        make
        
To compile and run the tests you issue the following command in the projectA directory:

        make tests
        
It is suggested not to use make tests as this project is still in development and there might be buggy functionalities in the test binary.
# Usage
All binaries are located in the bin directory.
To use them call them as follows:

        bin/projectA reference_graph.gfa node_list.txt
        bin/test
        
The reference_graph.gfa file should containe a reference graph in the .gfa format.
The node_list.txt file should containe clusterinformation appertaining to clusters on the reference graph. It should be formatted in the specific "node_list" format.
