#ifndef CSV
#define CSV

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>


std::vector<std::vector<std::string>> read_csv(const char* file){

    std::vector<std::vector<std::string>> df;

    std::ifstream inputFile(file); // Replace "example.txt" with the name of your file

    if (!inputFile.is_open()) {
        std::cerr << "Failed to open the file." << std::endl;
        return df;
    }

    std::string line;

    while (std::getline(inputFile, line)) {
        std::vector<std::string> row;
        // Process the line here
        // std::cout << line << std::endl;
        std::vector<std::string> tokens;
        std::istringstream tokenStream(line);
        std::string token;

        while (std::getline(tokenStream, token, ',')) {
            tokens.push_back(token);
        }

        // Now 'tokens' contains the individual parts of the line separated by commas
        for (int j = 0; j < tokens.size(); j++) {
            // std::cout << j << ": " << tokens[j] << std::endl;
            row.push_back(tokens[j]);
        }
        std::cout << std::endl;
        df.push_back(row);
    }

    inputFile.close(); // Close the file when you're done

    return df;

}

#endif