#include <fstream>
#include <iostream>
#include <vector>
#include <cstdio>

int main(int argc, char** argv)
{
    if(argc!=2){
        std::cerr << "Usage: read_packed filename" << std::endl;
        exit(1);
    }
    std::ifstream fin(argv[1], std::ios::in | std::ios::binary);
    std::vector<std::vector<unsigned short>> samples(128);
    uint32_t word=0;
    size_t nwords=0;
    uint32_t words[3]; // The shuffle pattern repeats every three words
    size_t index=0; // The index within words for the next word
    do{

        fin.read((char*)&word, sizeof(uint32_t));
        ++nwords;
        if(word==0xdeadbeef){
            // We got start-of-frame. Next word is the tdc
            fin.read((char*)&word, sizeof(uint32_t));
            std::cout << "Got new frame with TDC " << word << std::endl;
            if(index!=0){
                std::cout << "Error: index at new frame is " << index << ", not 0" << std::endl;
            }
            ++nwords;
            continue;
        }
        words[index]=word;
        if(index==2){
            // We've got our third word. Do the unpacking
            unsigned short thisSample[8];
            thisSample[0]=words[0] & 0xfff;
            thisSample[1]=(words[0] >> 12) & 0xfff;
            thisSample[2]=((words[0] >> 24) & 0xff) | ((words[1] << 8) & 0xf00);
            thisSample[3]=(words[1] >> 4) & 0xfff;
            thisSample[4]=(words[1] >> 16) & 0xfff;
            thisSample[5]=((words[1] >> 28) & 0xf) | ((words[2] << 4) & 0xff0);
            thisSample[6]=(words[2] >> 8) & 0xfff;
            thisSample[7]=(words[2] >> 20) & 0xfff;

            printf("% 6d % 6d % 6d % 6d % 6d % 6d % 6d % 6d\n",
                   thisSample[0], thisSample[1],
                   thisSample[2], thisSample[3],
                   thisSample[4], thisSample[5], 
                   thisSample[6], thisSample[7]);
            
            // Reset the index for the next set
            index=0;
        }
        else{
            ++index;
        }
    }
    while(word!=0xffffffff && nwords<1000);
    std::cout << "Exiting after reading " << nwords << " words" << std::endl;
}
