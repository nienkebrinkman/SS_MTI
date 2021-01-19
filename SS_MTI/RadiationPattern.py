import numpy as np


def R_P(take_off_angle, strike, dip, rake, az):
    """ Radiation pattern for P"""
    inc = np.deg2rad(take_off_angle)
    SR = Fault_geom_SR(dip, rake)
    QR = Fault_geom_QR(strike, dip, rake, az)
    PR = Fault_geom_PR(strike, dip, rake, az)

    RP = SR * (3 * np.cos(inc) ** 2 - 1) - QR * np.sin(2 * inc) - PR * np.sin(inc) ** 2
    return RP


def R_SV(take_off_angle, strike, dip, rake, az):
    """ Radiation pattern for SV"""
    inc = np.deg2rad(take_off_angle)
    SR = Fault_geom_SR(dip, rake)
    QR = Fault_geom_QR(strike, dip, rake, az)
    PR = Fault_geom_PR(strike, dip, rake, az)

    RSV = (3 / 2) * SR * np.sin(2 * inc) + QR * np.cos(2 * inc) + (1 / 2) * PR * np.sin(2 * inc)
    return RSV


def R_SH(take_off_angle, strike, dip, rake, az):
    """ Radiation pattern for SH"""
    inc = np.deg2rad(take_off_angle)
    QL = Fault_geom_QL(strike, dip, rake, az)
    PL = Fault_geom_PL(strike, dip, rake, az)

    RSH = -QL * np.cos(inc) - PL * np.sin(inc)
    return RSH


def Fault_geom_SR(dip, rake):
    """  Fault geometry factor for P - SV waves"""
    delta = np.deg2rad(dip)
    lambd = np.deg2rad(rake)

    SR = np.sin(lambd) * np.sin(delta) * np.cos(delta)
    return SR


def Fault_geom_QR(strike, dip, rake, az):
    """  Fault geometry factor for P - SV waves"""
    delta = np.deg2rad(dip)
    lambd = np.deg2rad(rake)

    phi_az = np.deg2rad(strike - az)

    QR = np.sin(lambd) * np.cos(2 * delta) * np.sin(phi_az) + np.cos(lambd) * np.cos(
        delta
    ) * np.cos(phi_az)
    return QR


def Fault_geom_PR(strike, dip, rake, az):
    """  Fault geometry factor for P - SV waves"""
    delta = np.deg2rad(dip)
    lambd = np.deg2rad(rake)

    phi_az = np.deg2rad(strike - az)

    PR = np.cos(lambd) * np.sin(delta) * np.sin(2 * phi_az) - np.sin(lambd) * np.sin(
        delta
    ) * np.cos(delta) * np.cos(2 * phi_az)
    return PR


def Fault_geom_PL(strike, dip, rake, az):
    """  Fault geometry factor for SH waves"""
    delta = np.deg2rad(dip)
    lambd = np.deg2rad(rake)

    phi_az = np.deg2rad(strike - az)

    PL = np.sin(lambd) * np.sin(delta) * np.cos(delta) * np.sin(2 * phi_az) + np.cos(
        lambd
    ) * np.sin(delta) * np.cos(2 * phi_az)
    return PL


def Fault_geom_QL(strike, dip, rake, az):
    """  Fault geometry factor for SH waves"""
    delta = np.deg2rad(dip)
    lambd = np.deg2rad(rake)

    phi_az = np.deg2rad(strike - az)

    QL = -np.cos(lambd) * np.cos(delta) * np.sin(phi_az) + np.sin(lambd) * np.cos(
        2 * delta
    ) * np.cos(phi_az)
    return QL
