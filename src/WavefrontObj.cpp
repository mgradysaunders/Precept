#include <iostream>
#include <fstream>
#include <sstream>
#include <pre-graphics/WavefrontObj>

namespace pre {

void WavefrontObj::read(const std::string& filename) {
    std::ifstream ifs(filename);
    if (!ifs.is_open())
        throw std::runtime_error("in WavefrontObj::read(), can't open file");
    read(ifs);
}

void WavefrontObj::read(std::istream& stream) {
    clear();
    std::string line;
    std::stringstream liness;
    long lineno = 1;
    auto next = [&]() {
        if (std::getline(stream, line)) {
            while (!line.empty() && std::isspace(line.back()))
                line.pop_back();
            while (!line.empty() && std::isspace(line.front()))
                line.erase(line.begin());
            lineno++;
            return true;
        }
        else
            return false;
    };
    auto starts_with = [&](std::string_view what) {
        if (line.starts_with(what)) {
            liness = std::stringstream(line.substr(what.size()));
            return true;
        }
        else
            return false;
    };
    std::uint16_t material = 0;
    std::map<std::string, std::uint16_t> material_indexes;
    material_indexes["default"] = 0;
    material_names[0] = "default";
    while (next()) {
        const char* failure = nullptr;
        if (starts_with("#"))
            continue; // Skip comments
        else if (starts_with("v ")) {
            auto& v = positions.v.emplace_back();
            liness >> v[0];
            liness >> v[1];
            liness >> v[2];
            if (not liness)
                failure = "failed to read v coordinate";
        }
        else if (starts_with("vt ")) {
            auto& v = texcoords.v.emplace_back();
            liness >> v[0];
            liness >> v[1];
            if (not liness)
                failure = "failed to read vt coordinate";
        }
        else if (starts_with("vn ")) {
            auto& v = normals.v.emplace_back();
            liness >> v[0];
            liness >> v[1];
            liness >> v[2];
            if (not liness)
                failure = "failed to read vn coordinate";
        }
        else if (starts_with("f ")) {
            std::string face;
            std::stringstream facess;
            auto maybe_read_then_consume = [&](int n) {
                int value = 0;
                char dummy;
                (facess >> value).clear();
                (facess >> dummy).clear();
                if (value < 0)
                    return value + n;
                else
                    return value - 1;
            };
            int sz = 0;
            while (liness >> face) {
                facess = std::stringstream(face);
                int fp = maybe_read_then_consume(positions.v.size());
                int ft = maybe_read_then_consume(texcoords.v.size());
                int fn = maybe_read_then_consume(normals.v.size());
                if (fp < 0)
                    failure = "failed to read v index";
                else if (ft < 0 and texcoords.f.size() > 0)
                    failure = "failed to read vt index (must be present for "
                              "every face or not at all)";
                else if (fn < 0 and normals.f.size() > 0)
                    failure = "failed to read vn index (must be present for "
                              "every face or not at all)";
                positions.f.emplace_back(fp);
                if (ft >= 0)
                    texcoords.f.emplace_back(ft);
                if (fn >= 0)
                    normals.f.emplace_back(fn);
                sz++;
            }
            if (sz < 3)
                failure = "face has less than 3 vertices";
            if (sz > 255)
                failure = "face has more than 255 vertices!";
            face_sizes.emplace_back(sz);
            face_materials.emplace_back(material);
        }
        else if (starts_with("usemtl ")) {
            material = material_indexes
                               .insert({liness.str(), material_indexes.size()})
                               .first->second;
            material_names[material] = liness.str();
        }
        if (failure) {
            clear();
            using namespace std::string_literals;
            throw std::runtime_error("in WavefrontObj::read(), "s + failure);
        }
    }
}

void WavefrontObj::write(const std::string& filename) const {
    std::ofstream ofs(filename);
    if (!ofs.is_open())
        throw std::runtime_error("in WavefrontObj::write(), can't open file");
    write(ofs);
}

void WavefrontObj::write(std::ostream& stream) const {
    int num_faces = face_sizes.size();
    if (num_faces == 0 or not positions)
        return;
    for (const auto& v : positions.v) {
        stream << "v ";
        stream << v[0] << ' ';
        stream << v[1] << ' ';
        stream << v[2] << '\n';
    }
    for (const auto& v : texcoords.v) {
        stream << "vt ";
        stream << v[0] << ' ';
        stream << v[1] << '\n';
    }
    for (const auto& v : normals.v) {
        stream << "vn ";
        stream << v[0] << ' ';
        stream << v[1] << ' ';
        stream << v[2] << '\n';
    }
    int offset = 0;
    for (int f = 0; f < num_faces; f++) {
        if (has_materials()) {
            if (f == 0 or face_materials[f - 1] != face_materials[f]) {
                stream << "usemtl ";
                auto itr = material_names.find(face_materials[f]);
                if (itr != material_names.end())
                    stream << itr->second;
                else
                    stream << face_materials[f];
                stream << '\n';
            }
        }
        stream << "f ";
        int num_verts = face_sizes[f];
        for (int v = 0; v < num_verts; v++) {
            stream << positions.f[offset] + 1;
            if (texcoords or normals) {
                stream << '/';
                if (texcoords)
                    stream << texcoords.f[offset] + 1;
            }
            if (normals) {
                stream << '/';
                stream << normals.f[offset] + 1;
            }
            stream << ' ';
            offset++;
        }
        stream << '\n';
    }
}

} // namespace pre
