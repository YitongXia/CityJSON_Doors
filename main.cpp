#include <iostream>
#include <fstream>
#include <string>
#include "include/json.hpp"
#include <cmath>

using json = nlohmann::json;

class Point{
public:
    double x;
    double y;
    double z;
    Point(double x1,double y1,double z1): x(x1),y(y1),z(z1)
    {}
};

// to select surfaces
static void parse_cityjson(json &j)
{
    for(auto &co:j["CityObjects"].items()) {
        for(auto &g:co.value()["geometry"]) {
            if(g["type"] == "Solid") {
                for(auto & surface:g["semantics"]["surfaces"]) {
                    std::cout<<surface<<std::endl;
                }
            }
        }
    }
}


// to select surfaces
// in this case, the surface is
static void read_surface_list(json &j, int surface_id)
{
    for(auto &co:j["CityObjects"].items()) {
        for(auto &g:co.value()["geometry"]) {
            if(g["type"] == "Solid") {
                json surfaces_list = g["semantics"]["values"];
                for(auto &surface:surfaces_list){
                }
            }
        }
    }
}

void split_surface(json &j ){

    for (auto& co : j["CityObjects"].items()) {
        for (auto &g: co.value()["geometry"]) {

            if (g["type"] == "Solid") {
                int roof_index = 0;
                int surface_index = 0;
                for (int i = 0; i < g["semantics"]["surfaces"].size(); ++i) {
                    if (g["semantics"]["surfaces"][i]["type"] == "WallSurface")
                    {
                        roof_index = i;
                        surface_index++;
                    }
                    else surface_index++;
                }
                for (auto &values: g["semantics"]["values"]) {

                    for (int k = 0; k < values.size(); ++k) {
                        if (values[k] == roof_index) {
                            for (auto &shell: g["boundaries"]) {
                                g["semantics"]["surfaces"][surface_index]["type"] = "WallSurface";
                                g["semantics"]["surfaces"][surface_index]["area"] = 0.0;
                                g["semantics"]["surfaces"][surface_index]["BuildingPart_id"] = co.key();
                                values[k] = surface_index;
                                surface_index++;
                            }
                        }
                    }
                }
            }
        }
    }
}

int main(int argc, const char * argv[]) {


    std::ifstream input("../data/cityjson.json");
    json j;
    input >> j;
    input.close();

    //parse_cityjson(j);
    read_surface_list(j,1);



    return 0;
}
