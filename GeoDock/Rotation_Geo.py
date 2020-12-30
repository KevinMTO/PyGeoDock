# author : Kevin Mato
# date : 07/05/2020 15:30
# Location :
import numpy



#############################################################
# REMEMBER!! USE .TO_RAD() TO FOLLOW ORIGINAL CODE CONVENTION
#############################################################


def matrix(angle_rad, o_vector, target_vector):

    x_vector = target_vector[0]
    y_vector = target_vector[1]
    z_vector = target_vector[2]
    x_orig = o_vector[0]
    y_orig = o_vector[1]
    z_orig = o_vector[2]

    u = x_vector - x_orig
    v = y_vector - y_orig
    w = z_vector - z_orig

    u2 = u * u
    v2 = v * v
    w2 = w * w

    l2 = u * u + v * v + w * w

    l = numpy.sqrt(l2)

    sint = numpy.sin(angle_rad)
    cost = numpy.cos(angle_rad)

    one_minus_cost = 1 - numpy.cos(angle_rad)

    return numpy.array([
        [
            (u2 + (v2 + w2) * cost) / l2,
            (u * v * one_minus_cost - w * l * sint) / l2,
            (u * w * one_minus_cost + v * l * sint) / l2,
            ( (x_orig * (v2 + w2) - u * (y_orig * v + z_orig * w) ) * one_minus_cost +
                        (y_orig * w - z_orig * v) * l * sint )/l2
        ],

        [
            ( (u * v * one_minus_cost) + (w * l * sint) )/l2,
            (v2 + (u2 + w2) * cost )/l2,
            ( (v * w * one_minus_cost) - (u * l * sint) )/ l2,
            ( ( y_orig * (u2 + w2)  - v * (x_orig * u + z_orig * w) ) * one_minus_cost +
                        (z_orig * u - x_orig * w) * l * sint )/l2
        ],

        [
            ( u * w * one_minus_cost - v * l * sint) / l2,
            ( v * w * one_minus_cost + u * l * sint) / l2,
            ( w2 + (u2 + v2) * cost ) / l2,
            (  ( z_orig * (u2 + v2) - w * (x_orig * u + y_orig * v) )  * one_minus_cost +
                        (x_orig * v - y_orig * u) * l * sint) / l2
        ],

        [0, 0, 0, 1]
    ])


def rounding(poly, val=3):

    return numpy.round(poly, val)


############################################################################################################
