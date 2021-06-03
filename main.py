from sympy import diff, symbols
import math
x_1 = symbols('x1')
x_2 = symbols('x2')


x0 = [[-1.2, 0.0]]
h = 0.1
E = 0.000001

matrix_h = [[], []]
matrix_h_obr = [[], []]
matrix = []
res_matrix = []
lambdf = []
napravlenie = []


def fx_add(x1, x2):
    u = diff("(10*(x1-x2)**2+(x1-1)**2)**(1/4)", x_1, 1)
    return eval(str(u))

def fx_func(x1, x2):
    return (10*(x1-x2)**2+(x1-1)**2)**(1/4)


def left_1(x, y, h):
    return (fx_func(x, y) - fx_func(x - h, y)) / h


def right_1(x, y, h):
    return (fx_func(x + h, y) - fx_func(x, y)) / h


def center_1_x(x, y, h):
    return (fx_func(x + h, y) - fx_func(x - h, y)) / (2 * h)


def center_1_y(x, y, h):
    return (fx_func(x, y + h) - fx_func(x, y - h)) / (2 * h)


def center_2_xy(x, y, h):
    return (fx_func(x + h, y + h) - fx_func(x + h, y - h) - fx_func(x - h, y + h) + fx_func(x - h, y - h)) / (4 * h ** 2)


def center_2_x(x, y, h):
    return (fx_func(x + h, y) - fx_func(x, y) - fx_func(x, y) + fx_func(x - h, y)) / (h ** 2)


def center_2_y(x, y, h):
    return (fx_func(x, y + h) - fx_func(x, y) - fx_func(x, y) + fx_func(x, y - h)) / (h ** 2)


def det(matrix_h):
    matrix_h[0] = [round(center_2_x(x0[-1][0], x0[-1][1], h), 5), round(center_2_xy(x0[-1][0], x0[-1][1], h), 5)]
    matrix_h[1] = [round(center_2_xy(x0[-1][0], x0[-1][1], h), 5), round(center_2_y(x0[-1][0], x0[-1][1], h), 5)]
    return matrix_h[0][0] * matrix_h[1][1] - matrix_h[0][1] * matrix_h[1][0]


def obr(matrix_h, matrix_h_obr):
    detet = det(matrix_h)
    matrix_h[0][1] = -matrix_h[0][1]
    matrix_h[1][0] = -matrix_h[1][0]

    matrix_h[0][0], matrix_h[1][1] = matrix_h[1][1], matrix_h[0][0]

    matrix_h_obr[0] = [round(matrix_h[0][0] / detet, 5), round(matrix_h[0][1] / detet, 5)]
    matrix_h_obr[1] = [round(matrix_h[1][0] / detet, 5), round(matrix_h[1][1] / detet, 5)]
    return matrix_h_obr


def matrix_1(matrix):
    matrix = [round(center_1_x(x0[-1][0], x0[-1][1], h), 5), round(center_1_y(x0[-1][0], x0[-1][1], h), 5)]
    return matrix


def mnoj(matrix_h, matrix, res_matrix, matrix_h_obr):
    matrix_h_obr = obr(matrix_h, matrix_h_obr)
    matrix = matrix_1(matrix)
    res_matrix = [round(matrix_h_obr[0][0] * matrix[0] + matrix_h_obr[0][1] * matrix[1], 5), round(matrix_h_obr[1][0] * matrix[0] + matrix_h_obr[1][1] * matrix[1], 5)]
    return res_matrix


def norm(array):
    return math.sqrt(sum([value ** 2 for value in array]))


def l(x, s):
    return 0.001 * (norm(x) / norm(s))


def naprav():
    return mnoj(matrix_h, matrix, res_matrix, matrix_h_obr)


def get_new_x(x, l, s):
    return [x[0] + l * s[0], x[1] + l * s[1]]


def sven(x0, s):
    is_end = False
    lambd = l(x0, s)
    power = 1
    lambdf.append(lambd)
    values = []
    x_plus = get_new_x(x0, lambd, s)
    x_minus = get_new_x(x0, - lambd, s)
    a = fx_func(*x_plus)
    b = fx_func(*x_minus)
    c = fx_func(*x0)
    if fx_func(*x_plus) < fx_func(*x_minus):
        values.append({"l": lambd, "x": x_plus, "f(x)": fx_func(*x_plus)})
    else:
        values.append({"l": -lambd, "x": x_minus, "f(x)": fx_func(*x_minus)})

    while not is_end:
        new_l = values[-1]["l"] + 2 ** power * values[0]["l"]
        x = get_new_x(x0, new_l, s)
        fx = fx_func(*x)

        values.append({"l": new_l, "x": x, "f(x)": fx})
        power += 1

        if values[-2]["f(x)"] < values[-1]["f(x)"]:
            middle_l = (values[-1]["l"] + values[-2]["l"]) / 2
            middle_x = get_new_x(x0, middle_l, s)
            values.append({"l": middle_l, "x": middle_x, "f(x)": fx_func(*middle_x)})
            is_end = True

    if len(values) == 3:
        # return values[0], {"l": (values[1]["l"] + values[2]["l"])/2, "x": numpy.ndarray([(values[1]["x"][0] + values[2]["x"][0])/2, (values[1]["x"][1] + values[2]["x"][1])/2]), "f(x)": (values[1]["f(x)"] + values[2]["f(x)"])/2}
        return values[0], values[1]
    else:
        return values[-4], values[-1]


def dsk_powell(x0, interval, accuracy, s):
    def get_middle_value(inter):
        l = (inter[0]["l"] + inter[-1]["l"]) / 2
        x = get_new_x(x0, l, s),
        return {"l": l,
                "x": x[0],
                "f(x)": fx_func(*x[0])
                }

    def get_dsk_approximate_min(inter, s):
        middle_value = get_middle_value(inter)
        l_min = middle_value["l"] + (
                ((inter[0]["f(x)"] - inter[1]["f(x)"]) * abs(middle_value["l"] - inter[0]["l"])) /
                (2 * (inter[0]["f(x)"] - 2 * middle_value["f(x)"] + inter[1]["f(x)"])))

        x_min = get_new_x(inter[0]["x"], l_min - inter[0]["l"], s)
        fx_min = fx_func(*x_min)

        return {
            "l": l_min,
            "x": x_min,
            "f(x)": fx_min
        }

    def full_dsk_interval(inter):
        return {
            "left": inter[0],
            "middle": get_middle_value(inter),
            "approximate_min": get_dsk_approximate_min(inter, s),
            "right": inter[-1]
        }

    def get_powell_approximate_min(inter):
        a1 = (inter[1]["f(x)"] - inter[0]["f(x)"]) / (inter[1]["l"] - inter[0]["l"])
        a2 = (((inter[2]["f(x)"] - inter[0]["f(x)"]) / (inter[2]["l"] - inter[0]["l"])) - a1) / (
                inter[2]["l"] - inter[1]["l"])
        l_min = (inter[0]["l"] + inter[1]["l"]) / 2 - (a1 / (2 * a2))

        x_min = get_new_x(x0, l_min, s)
        fx_min = fx_func(*x_min)
        return {
            "l": l_min,
            "x": x_min,
            "f(x)": fx_min
        }

    def full_powell_interval(inter):
        sorted_interval = sorted(inter, key=lambda point: point["l"])
        return {
            "left": sorted_interval[0],
            "middle": sorted_interval[1],
            "approximate_min": get_powell_approximate_min(sorted_interval),
            "right": sorted_interval[2]
        }

    values = [full_dsk_interval(interval)]
    points = sorted([value for value in values[-1].values()], key=lambda value: value["l"])
    values.append(full_powell_interval([points[1], points[2], points[3]]))
    return values



napravka = [-naprav()[0], -naprav()[1]]
sven_interval_1 = sven(x0[-1], napravka)

precise_interval_1 = dsk_powell(x0[-1], sven_interval_1, E, napravka)

lamda = precise_interval_1[1]['approximate_min']["l"]

x0.append([x0[-1][0] - lamda * mnoj(matrix_h, matrix, res_matrix, matrix_h_obr)[0],
            x0[-1][1] - lamda * mnoj(matrix_h, matrix, res_matrix, matrix_h_obr)[1]])


while (abs(x0[-2][0] - x0[-1][0]) > 0.00001 and abs(x0[-2][1] - x0[-1][1]) > 0.00001):
    napravka = [-naprav()[0], -naprav()[1]]
    sven_interval_1 = sven(x0[-1], napravka)
    precise_interval_1 = dsk_powell(x0[-1], sven_interval_1, E, napravka)

    lamda = precise_interval_1[1]['approximate_min']["l"]

    x0.append([x0[-1][0] - lamda * mnoj(matrix_h, matrix, res_matrix, matrix_h_obr)[0],
               x0[-1][1] - lamda * mnoj(matrix_h, matrix, res_matrix, matrix_h_obr)[1]])

print(x0)
