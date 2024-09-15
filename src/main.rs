mod solver;

use glutin_window::GlutinWindow as Window;
use opengl_graphics::{GlGraphics, OpenGL};
use piston::event_loop::{EventSettings, Events};
use piston::input::{RenderArgs, RenderEvent, UpdateArgs, UpdateEvent};
use piston::window::WindowSettings;
use graphics::*;
use graphics::types::Color;
use piston::{Button, ButtonArgs, ButtonEvent, ButtonState, Key, MouseCursorEvent, PressEvent};
use rand::random;

use solver::{EulerSolver};

// Multiplier for line width and colour value
const DECAY_RATE: f32 = 0.95;
const ALPHA_CULL: f32 = 0.15;
const SCALE: f64 = 30.0;
const PLACE_RADIUS: f64 = 4.0;
const NUM_POINTS: usize = 5000;

enum DrawMode {
    Point,
    Line
}

const MODE: DrawMode = DrawMode::Point;

pub struct App {
    gl: GlGraphics, // OpenGL drawing backend.
    points: Vec<(f64, f64, [f32; 4])>, // Position and alpha, for time fading
    origin: (f64, f64),
    solver: EulerSolver<NUM_POINTS>,
}

impl App {
    fn render(&mut self, args: &RenderArgs) {
        self.origin = (args.window_size[0] / 2.0, args.window_size[1] / 2.0);
        self.clear(args, [0.0, 0.0, 0.0, 1.0]);

        match MODE {
            DrawMode::Point => {
                for i in 0..self.points.len() {
                    let p = self.points[i];
                    let coords = (p.0 * SCALE, p.1 * SCALE);
                    let col = p.2;
                    if !is_outside_viewport(coords.0, coords.1, &args.viewport()) {
                        self.draw_point(args, coords, col);
                    }
                }
                self.points.clear();
            },
            DrawMode::Line => {
                for i in NUM_POINTS..self.points.len() {
                    let p0 = self.points[i - NUM_POINTS];
                    let p1 = self.points[i];

                    let col1 = p0.2;
                    let col2 = p1.2;

                    let avg_col = [
                        (col1[0] + col2[0]) * 0.5,
                        (col1[1] + col2[1]) * 0.5,
                        (col1[2] + col2[2]) * 0.5,
                        (col1[3] + col2[3]) * 0.5,
                    ];

                    if is_outside_viewport(p0.0, p0.1, &args.viewport()) || is_outside_viewport(p1.0, p1.1, &args.viewport()) {
                        continue;
                    }
                    else {
                        self.draw_line(args, (p0.0 * SCALE, p0.1 * SCALE), (p1.0 * SCALE, p1.1 * SCALE), avg_col);
                    }

                }
            }
        }
    }

    fn draw_point(&mut self, args: &RenderArgs, pos: (f64, f64), col: Color) {
        self.gl.draw(args.viewport(), |c, gl| {
            let circ = ellipse::circle(self.origin.0 + pos.0, self.origin.1 + pos.1, 0.5);
            ellipse(col, circ, c.transform, gl);
        });
    }

    fn draw_line(&mut self, args: &RenderArgs, pos1: (f64, f64), pos2: (f64, f64), col: [f32; 4]) {
        self.gl.draw(args.viewport(), |c, gl| {
            // let line = Line::new(col, col[3].into());

            let line = Line::new(col, (col[3]* 2.0).into());

            let x0 = self.origin.0 + pos1.0;
            let x1 = self.origin.0 + pos2.0;

            let y0 = self.origin.1 + pos1.1;
            let y1 = self.origin.1 + pos2.1;

            line.draw_from_to([x0, y0], [x1, y1], &c.draw_state, c.transform, gl);
        });
    }

    fn clear(&mut self, args: &RenderArgs, col: Color) {
        self.gl.draw(args.viewport(), |c, gl| {
            // Clear the screen.
            clear(col, gl);
        });
    }

    fn update(&mut self, args: &UpdateArgs) {
        self.solver.set_dt(args.dt * 10.0);
        self.solver.iterate();

        for (i, point) in self.solver.points().iter().enumerate() {
            let mut colour: [f32; 4] = [0.2, 1.0, 0.8, 1.0];

            self.points.push((point.0, point.1,
                              colour
            ));
        }

        match MODE {
            DrawMode::Point => {},
            DrawMode::Line => {
                let mut remove_list: Vec<usize> = Vec::new();

                // Fade points over time
                self.points.iter_mut().enumerate().for_each(|(i, p)| {
                    p.2[3] *= DECAY_RATE;
                    if p.2[3] <= ALPHA_CULL {
                        remove_list.push(i)
                    }
                    // p.2[2] *= 0.001 * DECAY_RATE;
                    // p.2[1] *= 0.05 * DECAY_RATE;

                });

                // Remove dead points
                remove_list.iter().rev().for_each(|i| {
                    self.points.remove(*i);
                });
            }
        }
    }

    fn handle_mouse(&mut self, x: f64, y: f64) {
        self.points.push((x - self.origin.0, y - self.origin.1, [1.0, 1.0, 1.0, 1.0]));
    }
}

fn main() {
    // Change this to OpenGL::V2_1 if not working.
    let opengl = OpenGL::V3_2;

    // Create a Glutin window.
    let mut window: Window = WindowSettings::new("parametric", [200, 200])
        .graphics_api(opengl)
        .exit_on_esc(true)
        .build()
        .unwrap();

    let mut x0 = [0.0; NUM_POINTS];
    let mut y0 = [0.0; NUM_POINTS];

    reset_points(&mut x0);
    reset_points(&mut y0);

    let mut solver = EulerSolver::<NUM_POINTS>::new(x0, y0, x_dot, y_dot, 0.0001);

    // Create a new game and run it.
    let mut app = App {
        gl: GlGraphics::new(opengl),
        points: Vec::new(),
        origin: (0.0, 0.0),
        solver: solver,
    };

    let mut events = Events::new(EventSettings::new());



    while let Some(e) = events.next(&mut window) {
        if let Some(args) = e.button_args() {
            if args.state == ButtonState::Press {
                match args.button {
                    Button::Keyboard(Key::R) => {
                        reset_points(&mut x0);
                        reset_points(&mut y0);
                        app.solver = EulerSolver::<NUM_POINTS>::new(x0, y0, x_dot, y_dot, 0.0001); // Unoptimised as of now
                    },
                    _ => ()
                }
            }
        }

        if let Some(args) = e.render_args() {
            app.render(&args);
        }

        if let Some(args) = e.update_args() {
            app.update(&args);
        }
    }
}

fn reset_points(x: &mut [f64; NUM_POINTS]) {
    for i in 0..NUM_POINTS {
        x[i] = (random::<f64>() * PLACE_RADIUS * 2.0) - PLACE_RADIUS;
    }
}

fn is_outside_viewport(x: f64, y: f64, viewport: &Viewport) -> bool {
    let [vx, vy, w, h] = viewport.rect;
    x < vx as f64 - w as f64 ||
    x > vx as f64 + w as f64 ||
    y < vy as f64 - h as f64 ||
    y > vy as f64 + h as f64
}

// Duffing
// const ALPHA: f64 = -1.0;
// const BETA: f64 = 1.0;
// const GAMMA: f64 = 4.0;
// const DELTA: f64 = 0.1;
// const OMEGA: f64 = 0.02;
//
// fn x_dot(x: f64, y: f64, t: f64) -> f64 {
//     GAMMA * (OMEGA * t).cos() - (DELTA * x) - (ALPHA * y) - (BETA * y * y * y)
// }
//
// fn y_dot(x: f64, y: f64, t: f64) -> f64 {
//     x
// }


// Van Der Pol
// const MU: f64 = 0.6;
//
// fn x_dot(x: f64, y: f64, t: f64) -> f64 {
//     MU * (x - 3.0_f64.recip() * x * x * x - y)
// }
//
// fn y_dot(x: f64, y: f64, t: f64) -> f64 {
//     MU.recip() * x
// }

// Thomas
// const ALPHA: f64 = 0.208186;
// fn x_dot(x: f64, y: f64, t: f64) -> f64 {
//     y.sin() - ALPHA * x
// }
//
// fn y_dot(x: f64, y: f64, t: f64) -> f64 {
//     x.sin() - ALPHA * y
// }

// Initial conditions, place radius: 3, Scale: 40, line mode


// Pendulum
const MU: f64 = 0.2;
const G: f64 = 0.981;
const L: f64 = 0.7;

fn x_dot(x: f64, y: f64, t: f64) -> f64 {
    y
}

fn y_dot(x: f64, y: f64, t: f64) -> f64 {
    -MU * y - (G / L) * x.sin()
}

// fn x_dot(x: f64, y: f64, t: f64) -> f64 {
//     y
// }
//
// fn y_dot(x: f64, y: f64, t: f64) -> f64 {
//     -0.2*y  - x*x*x.cos()*0.05
// }

// Initial conditions, place radius: 3, Scale: 60, line mode

// Halvorsen
// const ALPHA: f64 = 0.9;
// fn x_dot(x: f64, y: f64, t: f64) -> f64 {
//     - ALPHA * x - 4.0*y - y*y
// }
//
// fn y_dot(x: f64, y: f64, t: f64) -> f64 {
//     - ALPHA * y - 4.0*x - x*x
// }

// Sprott A
// fn x_dot(x: f64, y: f64, t: f64) -> f64 {
//    y
// }
//
// fn y_dot(x: f64, y: f64, t: f64) -> f64 {
//     -x + y*y -x*y
// }
// Use zoom 60, line, place radius 1

// FitzHugh-Nagumo
// const A: f64 = 0.2;
// const B: f64 = 0.03;
// const C: f64 = 2.0;
// fn x_dot(x: f64, y: f64, t: f64) -> f64 {
//     C * (y + x - x*x*x/3.0)
// }
//
// fn y_dot(x: f64, y: f64, t: f64) -> f64 {
//     -(x - A + B*y)/C
// }


// Duffing van der pol oscillator
// const A: f64 = 1.0; // Tweak for height of orbit
// const B: f64 = 0.03; // keep low
// const MU: f64 = 0.3; // Tweak for shape of attraction
// const G: f64 = 15.0; // Tweak for orbit band thickness
// const D: f64 = 0.01;
// fn x_dot(x: f64, y: f64, t: f64) -> f64 {
//     y
// }
//
// fn y_dot(x: f64, y: f64, t: f64) -> f64 {
//     -D*y - A*x -B*x*x*x + G*(0.1*t).cos() + MU*(1.0 - x*x)*y
// }

// fn x_dot(x: f64, y: f64, t: f64) -> f64 {
//     y
// }
//
// fn y_dot(x: f64, y: f64, t: f64) -> f64 {
//     -0.3 * y - 2.35 * x.sin() * (0.0005*t).sin()*x*y - 0.00005*t*x + 0.01*(1.0 - x*x)*y
// }
