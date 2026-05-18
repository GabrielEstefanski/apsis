//! Round-trip resume_state ↔ restore_resume_state for the three
//! integrators that publish non-empty payloads.

use crate::core::system::System;
use crate::domain::body::Body;
use crate::physics::integrator::IntegratorKind;
use crate::physics::integrator::ias15::Ias15;
use crate::physics::integrator::mercurius::Mercurius;
use crate::physics::integrator::traits::Integrator;
use crate::physics::integrator::whfast::WHFast;
use crate::units::UnitSystem;

fn sun_earth() -> Vec<Body> {
    vec![
        Body::star(1.0).at(0.0, 0.0).with_velocity(0.0, 0.0),
        Body::rocky(3e-6).at(1.0, 0.0).with_velocity(0.0, 1.0),
    ]
}

fn warmed_system(kind: IntegratorKind) -> System {
    let mut sys =
        System::new(sun_earth(), UnitSystem::canonical()).with_integrator(kind).with_dt(1e-3);
    for _ in 0..50 {
        sys.step();
    }
    sys
}

#[test]
fn whfast_round_trip_yields_bit_equal_state() {
    let sys = warmed_system(IntegratorKind::WHFast);
    let mut original = WHFast::default();
    original.restore_resume_state(&sys.integrator.resume_state()).unwrap();
    let bytes_a = original.resume_state();

    let mut clone = WHFast::default();
    clone.restore_resume_state(&bytes_a).unwrap();
    let bytes_b = clone.resume_state();
    assert_eq!(bytes_a, bytes_b);
}

#[test]
fn ias15_round_trip_yields_bit_equal_state() {
    let sys = warmed_system(IntegratorKind::Ias15);
    let mut original = Ias15::new();
    original.restore_resume_state(&sys.integrator.resume_state()).unwrap();
    let bytes_a = original.resume_state();

    let mut clone = Ias15::new();
    clone.restore_resume_state(&bytes_a).unwrap();
    let bytes_b = clone.resume_state();
    assert_eq!(bytes_a, bytes_b);
}

#[test]
fn mercurius_round_trip_yields_bit_equal_state() {
    let sys = warmed_system(IntegratorKind::Mercurius);
    let mut original = Mercurius::new();
    original.restore_resume_state(&sys.integrator.resume_state()).unwrap();
    let bytes_a = original.resume_state();

    let mut clone = Mercurius::new();
    clone.restore_resume_state(&bytes_a).unwrap();
    let bytes_b = clone.resume_state();
    assert_eq!(bytes_a, bytes_b);
}

#[test]
fn truncated_payload_errors() {
    use crate::physics::integrator::traits::ResumeError;
    let sys = warmed_system(IntegratorKind::Ias15);
    let bytes = sys.integrator.resume_state();
    let mut clone = Ias15::new();
    let err = clone.restore_resume_state(&bytes[..bytes.len() - 1]).unwrap_err();
    assert_eq!(err, ResumeError::Truncated);
}

#[test]
fn wrong_magic_errors() {
    use crate::physics::integrator::traits::ResumeError;
    let mut clone = Ias15::new();
    let err = clone.restore_resume_state(b"X15\x01\x00\x00\x00\x00").unwrap_err();
    assert_eq!(err, ResumeError::UnsupportedFormat);
}
