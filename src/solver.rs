pub struct EulerSolver2 {
    x: f64,
    y: f64,
    dx: Box<dyn Fn(f64, f64, f64) -> f64>, // F(x, y, t)
    dy: Box<dyn Fn(f64, f64, f64) -> f64>,
    dt: f64,
    t: f64
}

impl EulerSolver2 {
    pub fn new(x: f64, y: f64, dx: impl Fn(f64, f64, f64) -> f64 + 'static, dy: impl Fn(f64, f64, f64) -> f64 + 'static, dt: f64) -> Self {
        Self {
            x, y,
            dx: Box::new(dx),
            dy: Box::new(dy),
            dt,
            t: 0.0
        }
    }

    pub fn point(&self) -> (f64, f64) {
        (self.x, self.y)
    }

    pub fn iterate(&mut self) {
        let new_x = self.x + (self.dt * (self.dx)(self.x, self.y, self.t));
        let new_y = self.y + (self.dt * (self.dy)(self.x, self.y, self.t));
        self.t += self.dt;
        self.x = new_x;
        self.y = new_y;
    }
}

pub struct EulerSolver<const N: usize> {
    x: [f64; N],
    y: [f64; N],
    dx: Box<dyn Fn(f64, f64, f64) -> f64>, // F(x, y, t)
    dy: Box<dyn Fn(f64, f64, f64) -> f64>,
    dt: f64,
    t: f64
}

impl<const N: usize> EulerSolver<N> {
    pub fn new(x: [f64; N], y: [f64; N], dx: impl Fn(f64, f64, f64) -> f64 + 'static, dy: impl Fn(f64, f64, f64) -> f64 + 'static, dt: f64) -> Self {
        Self {
            x, y,
            dx: Box::new(dx),
            dy: Box::new(dy),
            dt,
            t: 0.0
        }
    }

    pub fn points(&self) -> [(f64, f64); N] {
        let mut out: Vec<(f64, f64)> = vec![(0.0, 0.0); N];
        for i in 0..N {
            out[i] = (self.x[i], self.y[i])
        };
        out.try_into().expect("Couldn't collect to slice")
    }

    pub fn iterate(&mut self) {
        for i in 0..N {
            let new_x = self.x[i] + (self.dt * (self.dx)(self.x[i], self.y[i], self.t));
            let new_y = self.y[i] + (self.dt * (self.dy)(self.x[i], self.y[i], self.t));
            self.t += self.dt;
            self.x[i] = new_x;
            self.y[i] = new_y;
        }
    }

    pub fn set_dt(&mut self, dt: f64) {
        self.dt = dt;
    }
}