#include <eggx.h>
#include <vector>
#include <ios>
#include <cmath>

#include "vector.hpp"
#include "gui.hpp"

GUI::GUI(int window_length, int border, int wait_time, System sys) {

    this->window_length = window_length;
    this->border        = border;
    this->wait_time     = wait_time;

    double volume = pow(double(window_length), 3.0);
    double mass   = volume * sys.density;

    scale = pow(mass / double(sys.n_part), 1.0 / 3.0);

    win = gopen(
            window_length + border * 2,
            window_length + border * 2);

    layer(win, 0, 1);

    gsetbgcolor(win, bgdark.c_str());
}

void GUI::redraw(System sys) {
    gclr(win);

    newcolor(win, divider.c_str());

    drawrect(
            win,
            double(border),
            double(border),
            double(window_length),
            double(window_length));

    newcolor(win, primary.c_str());

    for(int i=0; i < sys.n_part; i++) {
        // Projection of the particle position onto a 2d plane
        fillcirc(
                win,
                sys.pos[i].coord[0] * scale + double(border),
                sys.pos[i].coord[1] * scale + double(border),
                5.0,
                5.0);
    }

    copylayer(win, 1, 0);

    if(wait_time)
        msleep(wait_time);
}

void GUI::close() {
    // Get keyboard input from user and close the window
    ggetch();
    gclose(win);
}
