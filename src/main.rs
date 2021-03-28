#[derive(Clone, Copy, Debug)]
struct FieldPoint {
    u: f64,
    v: f64,
    P: f64,
    E: f64,
    rho: f64,
    k: f64,
    Cp: f64,
}

struct Solution {
    field: Vec<FieldPoint>,
    t: f64,
}

const ATM: f64 = 101325.0;

struct InitConds {
    P: f64,
    T: f64,
    R: f64,
    k: f64,
    Cp: f64,
    rho: f64,
    E: f64,
}

impl InitConds {
    fn new(P: f64, T: f64, R: f64, k: f64) -> Self {
        let Cp = R * k / (k - 1.0);
        let rho = P / (R * T);
        let E = P / ((k - 1.0) * rho);
        InitConds{ P, T, R, k, Cp, rho, E }
    }
}

fn air_init_conds() -> InitConds {
    InitConds::new(ATM, 297.0, 287.0, 1.4)
}

fn fuel_init_conds() -> InitConds {
    InitConds::new(40.0 * ATM, 400.0, 310.0, 1.2)
}

fn main() {
    let air = air_init_conds();
    let fuel = fuel_init_conds();


    println!("Hello, world!");
}
