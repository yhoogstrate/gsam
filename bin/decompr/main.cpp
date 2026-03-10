// Including the necessary libraries.
#include <iostream>
#include <fstream>
#include <zlib.h>

// multi-member gzip
 
namespace FileDecompressor {
 
    /**
    * @class GzipDecompressor
    * Represents a utility for decompressing gzip files.
    */
    class GzipDecompressor {
    public:
        /**
        * Decompresses a gzip file.
        *
        * @param inputFilePath The path to the gzip file.
        * @param outputFilePath The path to save the decompressed file.
        * @return bool True if the decompression is successful, false otherwise.
        */
        static bool DecompressFile(const std::string& inputFilePath, const std::string& outputFilePath) {
            // Open the input gzip file in binary mode.
            gzFile inputFile = gzopen(inputFilePath.c_str(), "rb");
            if (inputFile == nullptr) {
                std::cerr << "Error opening input file: " << inputFilePath << std::endl;
                return false;
            }
 
            // Open the output file in binary mode.
            std::ofstream outputFile(outputFilePath, std::ios::binary);
            if (!outputFile) {
                std::cerr << "Error opening output file: " << outputFilePath << std::endl;
                gzclose(inputFile);
                return false;
            }
 
            // Decompress the gzip file and write the decompressed data to the output file.
            char buffer[1024];
            int bytesRead;
            while ((bytesRead = gzread(inputFile, buffer, sizeof(buffer))) > 0) {
                outputFile.write(buffer, bytesRead);
            }
            
            std::cout << "EOF: " << gzeof(inputFile) << std::endl;
 
            // Check if any error occurred during decompression.
            int error;
            const char* errorMessage = gzerror(inputFile, &error);
            if (error != Z_OK) {
                std::cerr << "Error decompressing file: " << errorMessage << std::endl;
                gzclose(inputFile);
                outputFile.close();
                std::remove(outputFilePath.c_str());
                return false;
            }
 
            // Close the input and output files.
            gzclose(inputFile);
            outputFile.close();
 
            return true;
        }
    };
}
 
int main() {
    // Example usage:
    std::string inputFilePath = "/tmp/fq/ABA1-2544_NGS19-J611_BHM33WDSXX_S63_L003_R1_001.fastq.gz";
    std::string outputFilePath = "data.tar";
 
    if (FileDecompressor::GzipDecompressor::DecompressFile(inputFilePath, outputFilePath)) {
        std::cout << "File decompressed successfully!" << std::endl;
    } else {
        std::cerr << "Failed to decompress file." << std::endl;
    }
 
    return 0;
}
