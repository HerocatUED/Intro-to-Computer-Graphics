#include "Labs/5-Visualization/tasks.h"

#include <numbers>
#include <string>
#include <algorithm>
#include <cstdlib>

using VCX::Labs::Common::ImageRGB;
namespace VCX::Labs::Visualization {

    bool cmp0(Car const & a, Car const & b) { return a.cylinders < b.cylinders; }
    bool cmp1(Car const & a, Car const & b) {return a.displacement < b.displacement;}
    bool cmp2(Car const & a, Car const & b) { return a.weight < b.weight; }
    bool cmp3(Car const & a, Car const & b) { return a.horsepower < b.horsepower; }
    bool cmp4(Car const & a, Car const & b) {return a.acceleration < b.acceleration; }
    bool cmp5(Car const & a, Car const & b) { return a.mileage < b.mileage; }
    bool cmp6(Car const & a, Car const & b) { return a.year < b.year; } 

    struct CoordinateStates {
        // your code here
        std::vector<Car> data;
        std::pair<Car, Car> range;
        Car                      delta;
        glm::vec4                color         = glm::vec4(0, 0.5, 0, 0.3);
        glm::vec4                line_color    = color;
        glm::vec4                outline_color = glm::vec4(0, 0, 0, 1);
        float                    width         = 1;
        float                    x[7]          = { 0.05, 0.20, 0.35, 0.50, 0.65, 0.80, 0.95 };
        float                    up            = 0.1;
        float                    down          = 0.95;
        int                      index         = 0;

        CoordinateStates(std::vector<Car> const & d):
            data(d) {
            range.first.cylinders    = 10000;
            range.first.displacement = 10000;
            range.first.weight       = 10000;
            range.first.horsepower   = 10000;
            range.first.acceleration = 10000;
            range.first.mileage      = 10000;
            range.first.year         = 10000;

            range.second.cylinders    = 0;
            range.second.displacement = 0;
            range.second.weight       = 0;
            range.second.horsepower   = 0;
            range.second.acceleration = 0;
            range.second.mileage      = 0;
            range.second.year         = 0;

            FindRange();
            SortData(0);
        }

        void FindRange() {
            for (int i = 0; i < data.size(); ++i) {
                range.first.cylinders    = std::min(range.first.cylinders, data[i].cylinders);
                range.first.displacement = std::min(range.first.displacement, data[i].displacement);
                range.first.weight       = std::min(range.first.weight, data[i].weight);
                range.first.horsepower   = std::min(range.first.horsepower, data[i].horsepower);
                range.first.acceleration = std::min(range.first.acceleration, data[i].acceleration);
                range.first.mileage      = std::min(range.first.mileage, data[i].mileage);
                range.first.year         = std::min(range.first.year, data[i].year);

                range.second.cylinders    = std::max(range.second.cylinders, data[i].cylinders);
                range.second.displacement = std::max(range.second.displacement, data[i].displacement);
                range.second.weight       = std::max(range.second.weight, data[i].weight);
                range.second.horsepower   = std::max(range.second.horsepower, data[i].horsepower);
                range.second.acceleration = std::max(range.second.acceleration, data[i].acceleration);
                range.second.mileage      = std::max(range.second.mileage, data[i].mileage);
                range.second.year         = std::max(range.second.year, data[i].year);
            }
            delta.cylinders    = range.second.cylinders - range.first.cylinders;
            delta.displacement = range.second.displacement - range.first.displacement;
            delta.weight       = range.second.weight - range.first.weight;
            delta.horsepower   = range.second.horsepower - range.first.horsepower;
            delta.acceleration = range.second.acceleration - range.first.acceleration;
            delta.mileage      = range.second.mileage - range.first.mileage;
            delta.year         = range.second.year - range.first.year;
        }

        std::vector<float> FindCoordinate(Car const & c) {
            std::vector<float> ans;
            ans.push_back(up + 0.05 + 0.75 * (range.second.cylinders - c.cylinders) / delta.cylinders);
            ans.push_back(up + 0.05 + 0.75 * (range.second.displacement - c.displacement) / delta.displacement);
            ans.push_back(up + 0.05 + 0.75 * (range.second.weight - c.weight) / delta.weight);
            ans.push_back(up + 0.05 + 0.75 * (range.second.horsepower - c.horsepower) / delta.horsepower);
            ans.push_back(up + 0.05 + 0.75 * (range.second.acceleration - c.acceleration) / delta.acceleration);
            ans.push_back(up + 0.05 + 0.75 * (range.second.mileage - c.mileage) / delta.mileage);
            ans.push_back(up + 0.05 + 0.75 * (range.second.year - c.year) / delta.year);
            return ans;
        }

        void SortData(int index) {
            switch (index) {
            case (0): sort(data.begin(), data.end(), cmp0); break;
            case (1): sort(data.begin(), data.end(), cmp1); break;
            case (2): sort(data.begin(), data.end(), cmp2); break;
            case (3): sort(data.begin(), data.end(), cmp3); break;
            case (4): sort(data.begin(), data.end(), cmp4); break;
            case (5): sort(data.begin(), data.end(), cmp5); break;
            case (6): sort(data.begin(), data.end(), cmp6); break;
            }
        }

        void Paint(Common::ImageRGB & input) {
            for (int i = 0; i < data.size(); ++i) {
                std::vector<float> y = FindCoordinate(data[i]);
                for (int j = 0; j < 6; ++j) {
                    DrawLine(input, line_color, glm::vec2(x[j], y[j]), glm::vec2(x[j + 1], y[j + 1]), width);
                }
                line_color.b += 0.001;
                line_color.g -= 0.001;
            }

            for (int i = 0; i < 7; ++i) {
                DrawLine(input, outline_color, glm::vec2(x[i], up), glm::vec2(x[i], down), 2 * width);
                if (i == index) DrawFilledRect(input, glm::vec4(0.8, 1, 0.8, 0.7), glm::vec2(x[i] - 0.01, up), glm::vec2(0.02, down - up));
                else DrawFilledRect(input, glm::vec4(0.8, 0.8, 0.8, 0.5), glm::vec2(x[i] - 0.01, up), glm::vec2(0.02, down - up));
                DrawRect(input, glm::vec4(1, 1, 1, 0.7), glm::vec2(x[i] - 0.01, up - 0.02), glm::vec2(0.02, down - up + 0.02), 2 * width);
            }

            PrintText(input, outline_color, glm::vec2(x[0], 0.03), 0.02, "cylinders");
            PrintText(input, outline_color, glm::vec2(x[1], 0.03), 0.02, "displacement");
            PrintText(input, outline_color, glm::vec2(x[2], 0.03), 0.02, "weight");
            PrintText(input, outline_color, glm::vec2(x[3], 0.03), 0.02, "horsepower");
            PrintText(input, outline_color, glm::vec2(x[4], 0.03), 0.02, "acceleration (0-60 mph)");
            PrintText(input, outline_color, glm::vec2(x[5], 0.03), 0.02, "mileage");
            PrintText(input, outline_color, glm::vec2(x[6], 0.03), 0.02, "year");

            PrintText(input, outline_color, glm::vec2(x[0], 0.08), 0.02, std::to_string((int) (range.second.cylinders + delta.cylinders / 15)));
            PrintText(input, outline_color, glm::vec2(x[1], 0.08), 0.02, std::to_string((int) (range.second.displacement + delta.displacement / 15)));
            PrintText(input, outline_color, glm::vec2(x[2], 0.08), 0.02, std::to_string((int) (range.second.weight + delta.weight / 15)));
            PrintText(input, outline_color, glm::vec2(x[3], 0.08), 0.02, std::to_string((int) (range.second.horsepower + delta.horsepower / 15)));
            PrintText(input, outline_color, glm::vec2(x[4], 0.08), 0.02, std::to_string((int) (range.second.acceleration + delta.acceleration / 15)));
            PrintText(input, outline_color, glm::vec2(x[5], 0.08), 0.02, std::to_string((int) (range.second.mileage + delta.mileage / 15)));
            PrintText(input, outline_color, glm::vec2(x[6], 0.08), 0.02, std::to_string((int) (range.second.year + delta.year / 15)));

            PrintText(input, outline_color, glm::vec2(x[0], 0.97), 0.02, std::to_string((int) (range.first.cylinders - delta.cylinders / 15)));
            PrintText(input, outline_color, glm::vec2(x[1], 0.97), 0.02, std::to_string((int) (range.first.displacement - delta.displacement / 15)));
            PrintText(input, outline_color, glm::vec2(x[2], 0.97), 0.02, std::to_string((int) (range.first.weight - delta.weight / 15)));
            PrintText(input, outline_color, glm::vec2(x[3], 0.97), 0.02, std::to_string((int) (range.first.horsepower - delta.horsepower / 15)));
            PrintText(input, outline_color, glm::vec2(x[4], 0.97), 0.02, std::to_string((int) (range.first.acceleration - delta.acceleration / 15)));
            PrintText(input, outline_color, glm::vec2(x[5], 0.97), 0.02, std::to_string((int) (range.first.mileage - delta.mileage / 15)));
            PrintText(input, outline_color, glm::vec2(x[6], 0.97), 0.02, std::to_string((int) (range.first.year - delta.year / 15)));
        }

        bool Update(InteractProxy const& proxy) {
            if (! proxy.IsHovering()) return false;
            if (proxy.IsClicking()) {  
                glm::vec2 pos = proxy.MousePos();
                if (pos.y < up || pos.y > down) return false;
                for (int i = 0; i < 7; ++i) {
                    if (pos.x > x[i]-0.01 && pos.x < x[i] + 0.01) {
                        index = i;
                        SortData(index);
                        line_color = color;
                        return true;
                    }
                }
            }
            return false;
        }
    };

    bool PaintParallelCoordinates(Common::ImageRGB & input, InteractProxy const & proxy, std::vector<Car> const & data, bool force) {
        // your code here
        static CoordinateStates states(data);    // initialize
        SetBackGround(input, glm::vec4(1));                   
        bool                    change = states.Update(proxy); // update according to user input
        if (! force && ! change) return false;                 // determine to skip repainting
        states.Paint(input);                                   // visualize
        return true;
    }

    void LIC(ImageRGB & output, Common::ImageRGB const & noise, VectorField2D const & field, int const & step) {
        // your code here
    }
}; // namespace VCX::Labs::Visualization