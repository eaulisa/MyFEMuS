#pragma once
#include <sstream>
#include <iomanip>
#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <cassert>
#include <cstring>
#include <cstdint>
#include <iostream>

// -------------------- Base64 encode --------------------
 auto base64_encode = [](const unsigned char* data, size_t len) {
        static const char* enc =
            "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
        std::string out;
        out.reserve(4 * ((len + 2) / 3));
        for (size_t i = 0; i < len;) {
            uint32_t a = i < len ? data[i++] : 0;
            uint32_t b = i < len ? data[i++] : 0;
            uint32_t c = i < len ? data[i++] : 0;
            uint32_t triple = (a << 16) | (b << 8) | c;
            out.push_back(enc[(triple >> 18) & 0x3F]);
            out.push_back(enc[(triple >> 12) & 0x3F]);
            out.push_back(enc[(triple >> 6) & 0x3F]);
            out.push_back(enc[triple & 0x3F]);
        }
        size_t mod = len % 3;
        if (mod) {
            out[out.size() - 1] = '=';
            if (mod == 1) out[out.size() - 2] = '=';
        }
        return out;
    };

// -------------------- Binary DataArray writer --------------------
template <typename T>
static void write_binary_array(std::ofstream &os,
                               const std::string &type,
                               const std::string &name,
                               int ncomp,
                               const std::vector<T> &data) {
    uint32_t nbytes = data.size() * sizeof(T);   // 32-bit length header
    std::vector<unsigned char> buffer(sizeof(uint32_t) + nbytes);
    std::memcpy(buffer.data(), &nbytes, sizeof(uint32_t));
    if (nbytes > 0)
        std::memcpy(buffer.data() + sizeof(uint32_t), data.data(), nbytes);

    std::string encoded = base64_encode(buffer.data(), buffer.size());

    os << "        <DataArray type=\"" << type << "\"";
    if (!name.empty()) os << " Name=\"" << name << "\"";
    os << " NumberOfComponents=\"" << ncomp << "\" format=\"binary\">\n";
    os << encoded << "\n";
    os << "        </DataArray>\n";
}


// -------------------- Debugging helper --------------------
static void check_sizes(const char* label, size_t nentries,
                        int ncomp, size_t expected) {
    std::cout << "[DEBUG] " << label
              << " entries=" << nentries
              << " ncomp=" << ncomp
              << " total=" << nentries * ncomp
              << " expected~" << expected
              << std::endl;
}

