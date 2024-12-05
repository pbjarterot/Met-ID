use crate::mass_match::{calculate_y, parse_equation, EquationCoefficients};

#[test]
fn testing_y_calculations() {
    let coeffs: EquationCoefficients = EquationCoefficients {
        a: 0.0,
        b: 0.0,
        c: 0.0,
    };
    let x: f64 = 300.0;

    assert_eq!(calculate_y(&coeffs, x), Some(300.0));

    let coeffs: EquationCoefficients = EquationCoefficients {
        a: 0.0,
        b: 0.0,
        c: 1.0,
    };
    let x: f64 = 300.0;

    assert_eq!(calculate_y(&coeffs, x), Some(299.9997));
}

#[test]
fn testing_parse_equation() {
    assert_eq!(
        parse_equation("0"),
        Some(EquationCoefficients {
            a: 0.0,
            b: 0.0,
            c: 0.0
        })
    );
    assert_eq!(
        parse_equation("1"),
        Some(EquationCoefficients {
            a: 0.0,
            b: 0.0,
            c: 1.0
        })
    );
    assert_eq!(parse_equation("Some"), None);
    assert_eq!(parse_equation("1;1,1"), None);
    assert_eq!(
        parse_equation("1,1,1"),
        Some(EquationCoefficients {
            a: 1.0,
            b: 1.0,
            c: 1.0
        })
    );
    assert_eq!(
        parse_equation("1, 1, 1"),
        Some(EquationCoefficients {
            a: 1.0,
            b: 1.0,
            c: 1.0
        })
    );
}
