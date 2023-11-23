#include <random>

#include <spdlog/spdlog.h>

#include "Labs/1-Drawing2D/tasks.h"

using VCX::Labs::Common::ImageRGB;

namespace VCX::Labs::Drawing2D {
    /******************* 1.Image Dithering *****************/
    void DitheringThreshold(
        ImageRGB &       output,
        ImageRGB const & input) {
        for (std::size_t x = 0; x < input.GetSizeX(); ++x)
            for (std::size_t y = 0; y < input.GetSizeY(); ++y) {
                glm::vec3 color = input[{ x, y }];
                output.SetAt({ x, y }, {
                                           color.r > 0.5 ? 1 : 0,
                                           color.g > 0.5 ? 1 : 0,
                                           color.b > 0.5 ? 1 : 0,
                                       });
            }
    }

    void DitheringRandomUniform(
        ImageRGB &       output,
        ImageRGB const & input) {
        // your code here:
        std::default_random_engine             random(time(NULL));
        std::uniform_real_distribution<double> dis(-0.5, 0.5);
        for (std::size_t x = 0; x < input.GetSizeX(); ++x)
            for (std::size_t y = 0; y < input.GetSizeY(); ++y) {
                double    d     = dis(random);
                glm::vec3 color = input[{ x, y }];
                output.SetAt({ x, y }, {
                                           color.r + d,
                                           color.g + d,
                                           color.b + d,
                                       });
            }
        DitheringThreshold(output, output);
    }

    void DitheringRandomBlueNoise(
        ImageRGB &       output,
        ImageRGB const & input,
        ImageRGB const & noise) {
        // your code here:
        for (std::size_t x = 0; x < input.GetSizeX(); ++x)
            for (std::size_t y = 0; y < input.GetSizeY(); ++y) {
                glm::vec3 color_noise = noise[{ x, y }];
                glm::vec3 color_input = input[{ x, y }];
                output.SetAt({ x, y }, {
                                           0.5 * (color_input.r + color_noise.r),
                                           0.5 * (color_input.g + color_noise.g),
                                           0.5 * (color_input.b + color_noise.b),
                                       });
            }
        DitheringThreshold(output, output);
    }

    void DitheringOrdered(
        ImageRGB &       output,
        ImageRGB const & input) {
        // your code here:
        for (std::size_t x = 0; x < input.GetSizeX(); ++x)
            for (std::size_t y = 0; y < input.GetSizeY(); ++y) {
                glm::vec3 color = input[{ x, y }];
                if (color.r >=0.9) {
                    output.SetAt({ 3 * x + 1, 3 * y + 1 }, {
                                                                1,
                                                                1,
                                                                1,
                                                           });
                    output.SetAt({ 3 * x , 3 * y  }, {
                                                                 1,
                                                                 1,
                                                                 1,
                                                           });
                    output.SetAt({ 3 * x , 3 * y + 1 }, {
                                                                 1,
                                                                 1,
                                                                 1,
                                                           });
                    output.SetAt({ 3 * x, 3 * y + 2 }, {
                                                                 1,
                                                                 1,
                                                                 1,
                                                           });
                    output.SetAt({ 3 * x + 1, 3 * y + 2 }, {
                                                                 1,
                                                                 1,
                                                                 1,
                                                           });
                    output.SetAt({ 3 * x + 2, 3 * y + 2 }, {
                                                                 1,
                                                                 1,
                                                                 1,
                                                           });
                    output.SetAt({ 3 * x + 2, 3 * y + 1 }, {
                                                                 1,
                                                                 1,
                                                                 1,
                                                           });
                    output.SetAt({ 3 * x + 2, 3 * y }, {
                                                                 1,
                                                                 1,
                                                                 1,
                                                           });
                    output.SetAt({ 3 * x + 1, 3 * y }, {
                                                                 1,
                                                                 1,
                                                                 1,
                                                           });
                } 
                else if (color.r >= 0.8 && color.r < 0.9) {
                    output.SetAt({ 3 * x + 1, 3 * y + 1 }, {
                                                                 0,
                                                                 0,
                                                                 0,
                                                           });
                    output.SetAt({ 3 * x, 3 * y }, {
                                                         1,
                                                         1,
                                                         1,
                                                   });
                    output.SetAt({ 3 * x, 3 * y + 1 }, {
                                                             1,
                                                             1,
                                                             1,
                                                       });
                    output.SetAt({ 3 * x, 3 * y + 2 }, {
                                                             1,
                                                             1,
                                                             1,
                                                       });
                    output.SetAt({ 3 * x + 1, 3 * y + 2 }, {
                                                                 1,
                                                                 1,
                                                                 1,
                                                           });
                    output.SetAt({ 3 * x + 2, 3 * y + 2 }, {
                                                                 1,
                                                                 1,
                                                                 1,
                                                           });
                    output.SetAt({ 3 * x + 2, 3 * y + 1 }, {
                                                                 1,
                                                                 1,
                                                                 1,
                                                           });
                    output.SetAt({ 3 * x + 2, 3 * y }, {
                                                             1,
                                                             1,
                                                             1,
                                                       });
                    output.SetAt({ 3 * x + 1, 3 * y }, {
                                                             1,
                                                             1,
                                                             1,
                                                       });
                } 
                else if (color.r >= 0.7 && color.r < 0.8) {
                    output.SetAt({ 3 * x + 1, 3 * y + 1 }, {
                                                                 0,
                                                                 0,
                                                                 0,
                                                           });
                    output.SetAt({ 3 * x, 3 * y }, {
                                                         1,
                                                         1,
                                                         1,
                                                   });
                    output.SetAt({ 3 * x, 3 * y + 1 }, {
                                                             1,
                                                             1,
                                                             1,
                                                       });
                    output.SetAt({ 3 * x, 3 * y + 2 }, {
                                                             1,
                                                             1,
                                                             1,
                                                       });
                    output.SetAt({ 3 * x + 1, 3 * y + 2 }, {
                                                                 1,
                                                                 1,
                                                                 1,
                                                           });
                    output.SetAt({ 3 * x + 2, 3 * y + 2 }, {
                                                                 1,
                                                                 1,
                                                                 1,
                                                           });
                    output.SetAt({ 3 * x + 2, 3 * y + 1 }, {
                                                                 1,
                                                                 1,
                                                                 1,
                                                           });
                    output.SetAt({ 3 * x + 2, 3 * y }, {
                                                             1,
                                                             1,
                                                             1,
                                                       });
                    output.SetAt({ 3 * x + 1, 3 * y }, {
                                                             0,
                                                             0,
                                                             0,
                                                       });
                } 
                else if (color.r >= 0.6 && color.r < 0.7) {
                    output.SetAt({ 3 * x + 1, 3 * y + 1 }, {
                                                                 0,
                                                                 0,
                                                                 0,
                                                           });
                    output.SetAt({ 3 * x, 3 * y }, {
                                                         1,
                                                         1,
                                                         1,
                                                   });
                    output.SetAt({ 3 * x, 3 * y + 1 }, {
                                                             1,
                                                             1,
                                                             1,
                                                       });
                    output.SetAt({ 3 * x, 3 * y + 2 }, {
                                                             1,
                                                             1,
                                                             1,
                                                       });
                    output.SetAt({ 3 * x + 1, 3 * y + 2 }, {
                                                                 1,
                                                                 1,
                                                                 1,
                                                           });
                    output.SetAt({ 3 * x + 2, 3 * y + 2 }, {
                                                                 1,
                                                                 1,
                                                                 1,
                                                           });
                    output.SetAt({ 3 * x + 2, 3 * y + 1 }, {
                                                                 0,
                                                                 0,
                                                                 0,
                                                           });
                    output.SetAt({ 3 * x + 2, 3 * y }, {
                                                             1,
                                                             1,
                                                             1,
                                                       });
                    output.SetAt({ 3 * x + 1, 3 * y }, {
                                                             0,
                                                             0,
                                                             0,
                                                       });
                } 
                else if (color.r >= 0.5 && color.r < 0.6) {
                    output.SetAt({ 3 * x + 1, 3 * y + 1 }, {
                                                                 0,
                                                                 0,
                                                                 0,
                                                           });
                    output.SetAt({ 3 * x, 3 * y }, {
                                                         1,
                                                         1,
                                                         1,
                                                   });
                    output.SetAt({ 3 * x, 3 * y + 1 }, {
                                                             1,
                                                             1,
                                                             1,
                                                       });
                    output.SetAt({ 3 * x, 3 * y + 2 }, {
                                                             1,
                                                             1,
                                                             1,
                                                       });
                    output.SetAt({ 3 * x + 1, 3 * y + 2 }, {
                                                                 0,
                                                                 0,
                                                                 0,
                                                           });
                    output.SetAt({ 3 * x + 2, 3 * y + 2 }, {
                                                                 1,
                                                                 1,
                                                                 1,
                                                           });
                    output.SetAt({ 3 * x + 2, 3 * y + 1 }, {
                                                                 0,
                                                                 0,
                                                                 0,
                                                           });
                    output.SetAt({ 3 * x + 2, 3 * y }, {
                                                             1,
                                                             1,
                                                             1,
                                                       });
                    output.SetAt({ 3 * x + 1, 3 * y }, {
                                                             0,
                                                             0,
                                                             0,
                                                       });
                } 
                else if (color.r >= 0.4 && color.r < 0.5) {
                    output.SetAt({ 3 * x + 1, 3 * y + 1 }, {
                                                                 0,
                                                                 0,
                                                                 0,
                                                           });
                    output.SetAt({ 3 * x, 3 * y }, {
                                                         1,
                                                         1,
                                                         1,
                                                   });
                    output.SetAt({ 3 * x, 3 * y + 1 }, {
                                                             1,
                                                             1,
                                                             1,
                                                       });
                    output.SetAt({ 3 * x, 3 * y + 2 }, {
                                                             0,
                                                             0,
                                                             0,
                                                       });
                    output.SetAt({ 3 * x + 1, 3 * y + 2 }, {
                                                                 0,
                                                                 0,
                                                                 0,
                                                           });
                    output.SetAt({ 3 * x + 2, 3 * y + 2 }, {
                                                                 1,
                                                                 1,
                                                                 1,
                                                           });
                    output.SetAt({ 3 * x + 2, 3 * y + 1 }, {
                                                                 0,
                                                                 0,
                                                                 0,
                                                           });
                    output.SetAt({ 3 * x + 2, 3 * y }, {
                                                             1,
                                                             1,
                                                             1,
                                                       });
                    output.SetAt({ 3 * x + 1, 3 * y }, {
                                                             0,
                                                             0,
                                                             0,
                                                       });
                } 
                else if (color.r >= 0.3 && color.r < 0.4) {
                    output.SetAt({ 3 * x + 1, 3 * y + 1 }, {
                                                                 0,
                                                                 0,
                                                                 0,
                                                           });
                    output.SetAt({ 3 * x, 3 * y }, {
                                                         1,
                                                         1,
                                                         1,
                                                   });
                    output.SetAt({ 3 * x, 3 * y + 1 }, {
                                                             1,
                                                             1,
                                                             1,
                                                       });
                    output.SetAt({ 3 * x, 3 * y + 2 }, {
                                                             0,
                                                             0,
                                                             0,
                                                       });
                    output.SetAt({ 3 * x + 1, 3 * y + 2 }, {
                                                                 0,
                                                                 0,
                                                                 0,
                                                           });
                    output.SetAt({ 3 * x + 2, 3 * y + 2 }, {
                                                                 1,
                                                                 1,
                                                                 1,
                                                           });
                    output.SetAt({ 3 * x + 2, 3 * y + 1 }, {
                                                                 0,
                                                                 0,
                                                                 0,
                                                           });
                    output.SetAt({ 3 * x + 2, 3 * y }, {
                                                             0,
                                                             0,
                                                             0,
                                                       });
                    output.SetAt({ 3 * x + 1, 3 * y }, {
                                                             0,
                                                             0,
                                                             0,
                                                       });
                } 
                else if (color.r >= 0.2 && color.r < 0.3) {
                    output.SetAt({ 3 * x + 1, 3 * y + 1 }, {
                                                                 0,
                                                                 0,
                                                                 0,
                                                           });
                    output.SetAt({ 3 * x, 3 * y }, {
                                                         0,
                                                         0,
                                                         0,
                                                   });
                    output.SetAt({ 3 * x, 3 * y + 1 }, {
                                                             1,
                                                             1,
                                                             1,
                                                       });
                    output.SetAt({ 3 * x, 3 * y + 2 }, {
                                                             0,
                                                             0,
                                                             0,
                                                       });
                    output.SetAt({ 3 * x + 1, 3 * y + 2 }, {
                                                                 0,
                                                                 0,
                                                                 0,
                                                           });
                    output.SetAt({ 3 * x + 2, 3 * y + 2 }, {
                                                                 1,
                                                                 1,
                                                                 1,
                                                           });
                    output.SetAt({ 3 * x + 2, 3 * y + 1 }, {
                                                                 0,
                                                                 0,
                                                                 0,
                                                           });
                    output.SetAt({ 3 * x + 2, 3 * y }, {
                                                             0,
                                                             0,
                                                             0,
                                                       });
                    output.SetAt({ 3 * x + 1, 3 * y }, {
                                                             0,
                                                             0,
                                                             0,
                                                       });
                } 
                else if (color.r >= 0.1 && color.r < 0.2) {
                    output.SetAt({ 3 * x + 1, 3 * y + 1 }, {
                                                                 0,
                                                                 0,
                                                                 0,
                                                           });
                    output.SetAt({ 3 * x, 3 * y }, {
                                                         0,
                                                         0,
                                                         0,
                                                   });
                    output.SetAt({ 3 * x, 3 * y + 1 }, {
                                                             1,
                                                             1,
                                                             1,
                                                       });
                    output.SetAt({ 3 * x, 3 * y + 2 }, {
                                                             0,
                                                             0,
                                                             0,
                                                       });
                    output.SetAt({ 3 * x + 1, 3 * y + 2 }, {
                                                                 0,
                                                                 0,
                                                                 0,
                                                           });
                    output.SetAt({ 3 * x + 2, 3 * y + 2 }, {
                                                                 0,
                                                                 0,
                                                                 0,
                                                           });
                    output.SetAt({ 3 * x + 2, 3 * y + 1 }, {
                                                                 0,
                                                                 0,
                                                                 0,
                                                           });
                    output.SetAt({ 3 * x + 2, 3 * y }, {
                                                             0,
                                                             0,
                                                             0,
                                                       });
                    output.SetAt({ 3 * x + 1, 3 * y }, {
                                                             0,
                                                             0,
                                                             0,
                                                       });
                } 
                else if (color.r <0.1) {
                    output.SetAt({ 3 * x + 1, 3 * y + 1 }, {
                                                                 0,
                                                                 0,
                                                                 0,
                                                           });
                    output.SetAt({ 3 * x, 3 * y }, {
                                                         0,
                                                         0,
                                                         0,
                                                   });
                    output.SetAt({ 3 * x, 3 * y + 1 }, {
                                                             0,
                                                             0,
                                                             0,
                                                       });
                    output.SetAt({ 3 * x, 3 * y + 2 }, {
                                                             0,
                                                             0,
                                                             0,
                                                       });
                    output.SetAt({ 3 * x + 1, 3 * y + 2 }, {
                                                                 0,
                                                                 0,
                                                                 0,
                                                           });
                    output.SetAt({ 3 * x + 2, 3 * y + 2 }, {
                                                                 0,
                                                                 0,
                                                                 0,
                                                           });
                    output.SetAt({ 3 * x + 2, 3 * y + 1 }, {
                                                                 0,
                                                                 0,
                                                                 0,
                                                           });
                    output.SetAt({ 3 * x + 2, 3 * y }, {
                                                             0,
                                                             0,
                                                             0,
                                                       });
                    output.SetAt({ 3 * x + 1, 3 * y }, {
                                                             0,
                                                             0,
                                                             0,
                                                       });
                }
            }
    }

    void DitheringErrorDiffuse(
        ImageRGB &       output,
        ImageRGB const & input) {
        // your code here:
        size_t dy = input.GetSizeY() + 2;
        // considering y-1 or y+1 may go beyond boundary
        //  errors[1] to errors[input.GetSizeY()] stands for y=0 to y=input.GetSizeY() - 1
        //  errors[0] stands for (x,0 - 1),errors[input.GetSizeY()+1] stans for (x,input.GetSizeY())
        float * errors_r = new float[2 * dy];
        float * errors_g = new float[2 * dy];
        float * errors_b = new float[2 * dy];
        for (int i = 0; i < 2 * dy; ++i) {
            errors_r[i] = 0;
            errors_g[i] = 0;
            errors_b[i] = 0;
        }
        int cnt = 0;
        float sum_r, sum_g, sum_b, error_r, error_g, error_b;
        for (std::size_t x = 0; x < input.GetSizeX(); ++x) {
            for (std::size_t y = 0; y < input.GetSizeY(); ++y) {
                glm::vec3 color = input[{ x, y }];
                sum_r = color.r + errors_r[y + 1 + cnt * dy];
                sum_g = color.g + errors_g[y + 1 + cnt * dy];
                sum_b = color.b + errors_b[y + 1 + cnt * dy];
                output.SetAt({ x, y }, {
                                           sum_r >= 0.5 ? 1 : 0,
                                           sum_g >= 0.5 ? 1 : 0,
                                           sum_b >= 0.5 ? 1 : 0,
                                       });
                error_r = sum_r >= 0.5 ? sum_r - 1 : sum_r;
                error_g = sum_g >= 0.5 ? sum_g - 1 : sum_r;
                error_b = sum_b >= 0.5 ? sum_b - 1 : sum_r;
                //(x,y+1)
                errors_r[y + 2 + cnt * dy] += error_r * 5 / 16;
                errors_g[y + 2 + cnt * dy] += error_g * 5 / 16;
                errors_b[y + 2 + cnt * dy] += error_b * 5 / 16;
                //(x+1,y-1)
                errors_r[y + (cnt + 1) % 2 * dy] += error_r * 3 / 16;
                errors_g[y + (cnt + 1) % 2 * dy] += error_g * 3 / 16;
                errors_b[y + (cnt + 1) % 2 * dy] += error_b * 3 / 16;
                //(x+1,y)
                errors_r[y + 1 + (cnt + 1) % 2 * dy] += error_r * 7 / 16;
                errors_g[y + 1 + (cnt + 1) % 2 * dy] += error_g * 7 / 16;
                errors_b[y + 1 + (cnt + 1) % 2 * dy] += error_b * 7 / 16;
                //(x+1,y+1)
                errors_r[y + 2 + (cnt + 1) % 2 * dy] += error_r * 1 / 16;
                errors_g[y + 2 + (cnt + 1) % 2 * dy] += error_g * 1 / 16;
                errors_b[y + 2 + (cnt + 1) % 2 * dy] += error_b * 1 / 16;
            }
            for (int i = cnt * dy; i < cnt * dy + dy; ++i) {
                errors_r[i] = 0;
                errors_g[i] = 0;
                errors_b[i] = 0;
            }
            ++cnt;
            cnt %= 2;
        }
        delete[] errors_r;
        delete[] errors_g;
        delete[] errors_b;
    }

    /******************* 2.Image Filtering *****************/
    void Blur(
        ImageRGB &       output,
        ImageRGB const & input) {
        // your code here:
        size_t max_x = input.GetSizeX(), max_y = input.GetSizeY();
        float  r, g, b;
        for (std::size_t x = 0; x < max_x; ++x)
            for (std::size_t y = 0; y < max_y; ++y) {
                glm::vec3 color = input[{ x, y }];
                r = color.r, g = color.g, b = color.b;
                if (x > 0) {
                    if (y > 0) {
                        color = input[{ x - 1, y - 1 }];
                        r += color.r;
                        g += color.g;
                        b += color.b;
                    }
                    color = input[{ x - 1, y }];
                    r += color.r;
                    g += color.g;
                    b += color.b;
                    if (y < max_y - 1) {
                        color = input[{ x - 1, y + 1 }];
                        r += color.r;
                        g += color.g;
                        b += color.b;
                    }
                }
                if (y > 0) {
                    color = input[{ x, y - 1 }];
                    r += color.r;
                    g += color.g;
                    b += color.b;
                }
                if (y < max_y - 1) {
                    color = input[{ x, y + 1 }];
                    r += color.r;
                    g += color.g;
                    b += color.b;
                }
                if (x < max_x - 1) {
                    if (y > 0) {
                        color = input[{ x + 1, y - 1 }];
                        r += color.r;
                        g += color.g;
                        b += color.b;
                    }
                    color = input[{ x + 1, y }];
                    r += color.r;
                    g += color.g;
                    b += color.b;
                    if (y < max_y - 1) {
                        color = input[{ x + 1, y + 1 }];
                        r += color.r;
                        g += color.g;
                        b += color.b;
                    }
                }
                r /= 9;
                g /= 9;
                b /= 9;
                output.SetAt({ x, y }, {
                                           r,
                                           g,
                                           b,
                                       });
            }
    }

    void Edge(
        ImageRGB &       output,
        ImageRGB const & input) {
        // your code here:
        size_t max_x = input.GetSizeX(), max_y = input.GetSizeY();
        float  r1, r2, g1, g2, b1, b2;
        glm::vec3 color;
        for (std::size_t x = 0; x < max_x; ++x)
            for (std::size_t y = 0; y < max_y; ++y) {
                r1 = 0;
                g1 = 0;
                b1 = 0;
                if (x > 0) {
                    if (y > 0) {
                        color = input[{ x - 1, y - 1 }];
                        r1 -= color.r;
                        g1 -= color.g;
                        b1 -= color.b;
                    }
                    color = input[{ x - 1, y }];
                    r1 -= 2 * color.r;
                    g1 -= 2 * color.g;
                    b1 -= 2 * color.b;
                    if (y < max_y - 1) {
                        color = input[{ x - 1, y + 1 }];
                        r1 -= color.r;
                        g1 -= color.g;
                        b1 -= color.b;
                    }
                }
                if (x < max_x - 1) {
                    if (y > 0) {
                        color = input[{ x + 1, y - 1 }];
                        r1 += color.r;
                        g1 += color.g;
                        b1 += color.b;
                    }
                    color = input[{ x + 1, y }];
                    r1 += 2 * color.r;
                    g1 += 2 * color.g;
                    b1 += 2 * color.b;
                    if (y < max_y - 1) {
                        color = input[{ x + 1, y + 1 }];
                        r1 += color.r;
                        g1 += color.g;
                        b1 += color.b;
                    }
                }
                r2 = 0;
                g2 = 0;
                b2 = 0;
                if (y > 0) {
                    if (x > 0) {
                        color = input[{ x - 1, y - 1 }];
                        r2 -= color.r;
                        g2 -= color.g;
                        b2 -= color.b;
                    }
                    color = input[{ x, y - 1 }];
                    r2 -= 2 * color.r;
                    g2 -= 2 * color.g;
                    b2 -= 2 * color.b;
                    if (x < max_x - 1) {
                        color = input[{ x + 1, y - 1 }];
                        r2 -= color.r;
                        g2 -= color.g;
                        b2 -= color.b;
                    }
                }
                if (y < max_y - 1) {
                    if (x > 0) {
                        color = input[{ x - 1, y + 1 }];
                        r2 += color.r;
                        g2 += color.g;
                        b2 += color.b;
                    }
                    color = input[{ x, y + 1 }];
                    r2 += 2 * color.r;
                    g2 += 2 * color.g;
                    b2 += 2 * color.b;
                    if (x < max_x - 1) {
                        color = input[{ x + 1, y + 1 }];
                        r2 += color.r;
                        g2 += color.g;
                        b2 += color.b;
                    }
                }
                r1 = r1 > 0 ? r1 : -r1;
                r2 = r2 > 0 ? r2 : -r2;
                g1 = g1 > 0 ? g1 : -g1;
                g2 = g2 > 0 ? g2 : -g2;
                b1 = b1 > 0 ? b1 : -b1;
                b2 = b2 > 0 ? b2 : -b2;
                output.SetAt({ x, y }, {
                                           r1 > r2 ? r1 : r2,
                                           g1 > g2 ? g1 : g2,
                                           b1 > b2 ? b1 : b2,
                                       });
            }
    }

    /******************* 3. Image Inpainting *****************/
    void Inpainting(
        ImageRGB &         output,
        ImageRGB const &   inputBack,
        ImageRGB const &   inputFront,
        const glm::ivec2 & offset) {
        output             = inputBack;
        size_t      width  = inputFront.GetSizeX();
        size_t      height = inputFront.GetSizeY();
        glm::vec3 * g      = new glm::vec3[width * height];
        memset(g, 0, sizeof(glm::vec3) * width * height);
        // set boundary condition
        for (std::size_t y = 0; y < height; ++y) {
            // set boundary for (0, y), your code: g[y * width] = ?
            g[y * width] = inputBack[{ (size_t) offset.x, (size_t) offset.y + y }] - inputFront[{ 0, y }];
            // set boundary for (width - 1, y), your code: g[y * width + width - 1] = ?
            g[y * width + width - 1] = inputBack[{ width - 1 + (size_t) offset.x, y + (size_t) offset.y }] - inputFront[{ width - 1, y }];
        }
        for (std::size_t x = 0; x < width; ++x) {
            // set boundary for (x, 0), your code: g[x] = ?
            g[x] = inputBack[{ (size_t) offset.x + x, (size_t) offset.y }] - inputFront[{ x, 0 }];
            // set boundary for (x, height - 1), your code: g[(height - 1) * width + x] = ?
            g[(height - 1) * width + x] = inputBack[{ (size_t) offset.x + x, (size_t) offset.y + height - 1 }] - inputFront[{ x, height - 1 }];
        }

        // Jacobi iteration, solve Ag = b
        for (int iter = 0; iter < 8000; ++iter) {
            for (std::size_t y = 1; y < height - 1; ++y)
                for (std::size_t x = 1; x < width - 1; ++x) {
                    g[y * width + x] = (g[(y - 1) * width + x] + g[(y + 1) * width + x] + g[y * width + x - 1] + g[y * width + x + 1]);
                    g[y * width + x] = g[y * width + x] * glm::vec3(0.25);
                }
        }

        for (std::size_t y = 0; y < inputFront.GetSizeY(); ++y)
            for (std::size_t x = 0; x < inputFront.GetSizeX(); ++x) {
                glm::vec3 color = g[y * width + x] + inputFront.GetAt({ x, y });
                output.SetAt({ x + offset.x, y + offset.y }, color);
            }
        delete[] g;
    }

    /******************* 4. Line Drawing *****************/
    void DrawLine(
        ImageRGB &       canvas,
        glm::vec3 const  color,
        glm::ivec2 const p0,
        glm::ivec2 const p1) {
        // your code here:
        glm::ivec2 p2    = p0.x <= p1.x ? p0 : p1;
        glm::ivec2 p3    = p0.x > p1.x ? p0 : p1;
        float      slope = (float) (p2.y - p3.y) / (p2.x - p3.x);
        if (slope >= 0 && slope <= 1) {
            int x = p2.x, y = p2.y;
            int dx = 2 * (p3.x - p2.x), dy = 2 * (p3.y - p2.y);
            int dydx = dy - dx, F = dy - dx / 2;
            for (x; x <= p3.x; ++x) {
                canvas.SetAt({ (size_t) x, (size_t) y }, color);
                if (F < 0) F += dy;
                else {
                    ++y;
                    F += dydx;
                }
            }
        } else if (slope > 1) {
            int x = p2.x, y = p2.y;
            int dx = 2 * (p3.x - p2.x), dy = 2 * (p3.y - p2.y);
            int dydx = dx - dy, F = dx - dy / 2;
            for (y; y <= p3.y; ++y) {
                canvas.SetAt({ (size_t) x, (size_t) y }, color);
                if (F < 0) F += dx;
                else {
                    ++x;
                    F += dydx;
                }
            }
        } else if (slope < 0 && slope >= -1) {
            int x = p2.x, y = p2.y;
            int dx = 2 * (p3.x - p2.x), dy = 2 * (p2.y - p3.y);
            int dydx = dy - dx, F = dy - dx / 2;
            for (x; x <= p3.x; ++x) {
                canvas.SetAt({ (size_t) x, (size_t) y }, color);
                if (F < 0) F += dy;
                else {
                    --y;
                    F += dydx;
                }
            }
        } else if (slope < -1) {
            int x = p2.x, y = p2.y;
            int dx = 2 * (p3.x - p2.x), dy = 2 * (p2.y - p3.y);
            int dydx = dx - dy, F = dx - dy / 2;
            for (y; y >= p3.y; --y) {
                canvas.SetAt({ (size_t) x, (size_t) y }, color);
                if (F < 0) F += dx;
                else {
                    ++x;
                    F += dydx;
                }
            }
        }
    }

    /******************* 5. Triangle Drawing *****************/
    void DrawTriangleFilled(
        ImageRGB &       canvas,
        glm::vec3 const  color,
        glm::ivec2 const p0,
        glm::ivec2 const p1,
        glm::ivec2 const p2) {
        // your code here:
        glm::ivec2 p3, p4, p5;
        // sort points, p3.y >= p4.y >= p5.y
        p3 = p0.y > p1.y ? p0 : p1;
        p3 = p3.y > p2.y ? p3 : p2;
        p5 = p0.y <= p1.y ? p0 : p1;
        p5 = p5.y <= p2.y ? p5 : p2;
        if (p0.x != p3.x && p0.x != p5.x && p0.y != p3.y && p0.y != p5.y) p4 = p0;
        else if (p1.x != p3.x && p1.x != p5.x && p1.y != p3.y && p1.y != p5.y) p4 = p1;
        else p4 = p2;
        size_t  y  = p5.y, x;
        double  x0 = p5.x, x1 = p5.x;
        double  slope0 = (x0 - p3.x) / ((double) y - p3.y);
        double  slope1 = (x1 - p4.x) / ((double) y - p4.y);
        double &xL = x0 + slope0 < x1 + slope1 ? x0 : x1, &xR = x0 + slope0 > x1 + slope1 ? x0 : x1;
        for (y; y < p4.y; ++y) {
            for (x = (size_t) xL; x <= xR; ++x) {
                canvas.SetAt({ x, y }, color);
            }
            x0 += slope0;
            x1 += slope1;
        }
        slope0 = (x0 - p3.x) / ((double) y - p3.y);
        slope1 = (x1 - p3.x) / ((double) y - p3.y);
        for (y; y <= p3.y; ++y) {
            for (x = (size_t) xL; x <= xR; ++x) {
                canvas.SetAt({ x, y }, color);
            }
            x0 += slope0;
            x1 += slope1;
        }
    }

    /******************* 6. Image Supersampling *****************/
    void Supersample(
        ImageRGB &       output,
        ImageRGB const & input,
        int              rate) {
        // your code here:
        double times_x = (double)input.GetSizeX() / output.GetSizeX(), times_y = (double)input.GetSizeY() / output.GetSizeY();
        for (std::size_t x = 0; x < output.GetSizeX(); ++x)
            for (std::size_t y = 0; y < output.GetSizeY(); ++y) {
                glm::vec3 color ;
                double    r = 0, g = 0, b = 0;
                for (int i = 0; i < rate; ++i) {
                    double tx = rand() % (int)times_x + x * times_x, ty = rand() % (int)times_y + y * times_y;
                    color = input[{ (size_t) tx, (size_t) ty }];
                    r += color.r;
                    g += color.g;
                    b += color.b;
                }
                output.SetAt({ x, y }, {
                                          r/rate ,
                                          g/rate ,
                                          b/rate ,
                                       });
            }
    }

    /******************* 7. Bezier Curve *****************/
    glm::vec2 CalculateBezierPoint(
        std::span<glm::vec2> points,
        float const          t) {
        // your code here:
        size_t len = points.size();
        float * xs = new float[len];
        float * ys  = new float[len];
        for (size_t i = 0; i < len; ++i) {
            xs[i] = points[i].x;
            ys[i] = points[i].y;
        }
        while (len != 1) {
            for (size_t i = 0; i < len-1; ++i) {
                xs[i] = (1 - t) * xs[i] + t * xs[i + 1];
                ys[i] = (1 - t) * ys[i] + t * ys[i + 1];
            }
            len--;
        }
        int point_x = xs[0], point_y = ys[0];
        glm::vec2 point = points.front();
        point.x         = point_x;
        point.y         = point_y;
        delete[] xs;
        delete[] ys;
        return point;
    }
} // namespace VCX::Labs::Drawing2D