#pragma once

#include <string>

#include "system.hpp"

class GUI {
    private:
        int win;
        int window_length;
        int border;
        int wait_time;

        double scale;

        std::string bgdark="#212121";
        std::string bglight="#B2DFDB";
        std::string accent="#00BCD4";
        std::string light="#BDBDBD";
        std::string primary="#009688";
        std::string divider="#BDBDBD";
        std::string darkprimary="#212121";

    public:
        GUI(int window_length, int border, int wait_time, System sys);

        void redraw(System sys);

        void close();
};
