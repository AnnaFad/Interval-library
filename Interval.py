import math
import numpy as np
import sympy
import decimal as dec

class Interval:
    def __init__(self, x):
        """
        Constructor for Interval

        Parameters
        ----------
        x - list with size 2 [a, b]. Where a < b. a - lower bound, b - upper bound.
        """
        if len(x) != 2:
            raise ValueError("Length of the list in parameters must be 2.")
        if x[0] > x[1]:
            return None
        self.x = x.copy()

    def __repr__(self):
        result = "["
        if isinstance(self.x[0], (float, sympy.Float, np.float32, np.float64)):
            result += str(round(float(self.x[0]), 3))
        else:
            result += str(self.x[0])
        result += ", "
        if isinstance(self.x[1], (float, sympy.Float, np.float32, np.float64)):
            result += str(round(float(self.x[1]), 3))
        else:
            result += str(self.x[1])
        result += "]"
        return result

    def mid(self):
        """
        Returns central point of the interval
        """
        return 0.5 * (self.x[0] + self.x[1])


    def width(self):
        """
        Returns width of the interval
        """
        return self.x[1] - self.x[0]

    def radius(self):
        """
        Returns radius of the interval
        """
        return 0.5 * (self.x[0] - self.x[1])

    def scale(self, factor):
        """
        Changes width of the interval in given number of times with respect to the middle

        Parameters
        ----------
        factor - number of times in which interval's width is changing
        """
        m = self.mid()
        r = self.radius()
        self.x[0] = m - factor * r
        self.x[1] = m + factor * r

    def is_in(self, other):
        """
        Checks, is second interval or value in this interval.

        Parameters
        ----------
        other - interval or value

        Returns:
        -------
        True if other lies in self, False otherwise
        """

        ointerval = value_to_interval(other)
        return (self.x[0] >= ointerval.x[0]) and (self.x[1] <= ointerval.x[1])

    def is_no_intersec(self, other):
        """
        Checks that two intervals has no intersection

        Parameters
        ----------
        other - interval or value

        Returns:
        -------
        True if there is no intersetion of the intervals, False otherwise
        """
        ointerval = value_to_interval(other)
        return (self.x[0] > ointerval.x[1]) or (self.x[1] < ointerval.x[0])

    def intersec(self, other):
        """
        Returns intersection of two intervals

        Parameters
        ----------
        other - interval or value
        """
        ointerval = value_to_interval(other)
        if self.x[1] < ointerval.x[0] or ointerval.x[1] < self.x[0]:
            return None
        return Interval([max(self.x[0], ointerval.x[0]), min(self.x[1], ointerval.x[1])])

    def __getitem__(self, item):
        return self.x[item]

    def __setitem__(self, key, value):
        self.x.__setitem__(key, value)

    def __neg__(self):
        ninterval = Interval(self.x)
        ninterval.x[0] = - self.x[1]
        ninterval.x[1] = - self.x[0]
        return ninterval

    def __add__(self, other):
        ointerval = value_to_interval(other)
        ninterval = Interval(self.x)
        ninterval.x[0] = self.x[0] + ointerval.x[0]
        ninterval.x[1] = self.x[1] + ointerval.x[1]
        return ninterval

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        ointerval = value_to_interval(other)
        return self.__add__(-ointerval)


    def __rsub__(self, other):
        ointerval = value_to_interval(other)
        return ointerval.__add__(-self)

    def __pow__(self, factor):
        if ((not isinstance(factor, int)) and (not isinstance(factor, sympy.Integer))) and ((not isinstance(factor, float)) and
         (not isinstance(factor, sympy.Float))):
            raise TypeError("Factor for exponentiation should be int or float")
        ninterval = Interval(self.x)
        is_nec_fact = False
        if factor < 0:
          factor *= -1
          is_nec_fact = True
        u = self.x[0] ** factor
        v = self.x[1] ** factor
        if factor == 0:
            ninterval.x[0] = 1
            ninterval.x[1] = 1
        elif factor % 2 == 0:
            ninterval.x[1] = max(u, v)
            if self.x[0] <= 0 and self.x[1] >= 0:
                ninterval.x[0] = 0
            else:
                ninterval.x[0] = min(u, v)
        else:
            ninterval.x[0] = min(u, v)
            ninterval.x[1] = max(u, v)
        if is_nec_fact:
            ninterval = 1/ninterval
        return ninterval

    def __rpow__(self, value):
        if ((not isinstance(value, int)) and (not isinstance(value, sympy.Integer))) and ((not isinstance(value, float)) and
         (not isinstance(value, sympy.Float))):
            raise TypeError("Value for exponentiation should be int or float")
        if value < 0:
            raise ValueError("Value should be positive or 0")
        if value == 0:
            if self.x[0] < 0:
                raise ValueError("0 can't be in negative power")
            elif self.x[0] == 0 and self.x[1] == 0:
                return Interval([1, 1])
            elif self.x[0] == 0 and self.x[1] > 0:
                return [Interval([1, 1]), Interval([0, 0])]
            else:
                u = value**self.x[0]
                v = value**self.x[1]
                return Interval([min(u, v), max(u, v)])


    def __mul__(self, other):
        ointerval = value_to_interval(other)
        v = [self.x[0] * ointerval.x[0], self.x[0] * ointerval.x[1], self.x[1] * ointerval.x[0], self.x[1] * ointerval.x[1]]
        b = [min(v), max(v)]
        return Interval(b)


    def __my_true_div(self, a, b):
        if b == float("inf") or b == float("-inf") or b == dec.Decimal('Infinity') or b == dec.Decimal('-Infinity'):
            if a == float("inf") or a == float("-inf") or a == dec.Decimal('Infinity') or a == dec.Decimal('-Infinity'):
                raise ValueError('Infinity by infinity division')
            else:
                return 0
        return a / b

    def __my_floor_div(self, a, b):
        if b == float("inf") or b == float("-inf") or b == dec.Decimal('Infinity') or b == dec.Decimal('-Infinity'):
            if a == float("inf") or a == float("-inf") or a == dec.Decimal('Infinity') or a == dec.Decimal('-Infinity'):
                raise ValueError('Infinity by infinity division')
            else:
                return 0
        return a // b

    def __truediv__(self, other):
        ointerval = value_to_interval(other)
        if ointerval.x[0] == ointerval.x[1] == 0:
            if self.x[0] <= 0 <= self.x[1]:
                return Interval([float('-inf'), float('inf')])
            else:
                return None
        elif ointerval.x[0] > 0 or ointerval.x[1] < 0:
            return self.__mul__(Interval([1 / ointerval.x[1], 1 / ointerval.x[0]]))
        elif self.x[0] <= 0 <= self.x[1]:
            return Interval([float('-inf'), float('inf')])
        elif ointerval.x[0] == 0:
            if self.x[1] < 0:
                return Interval([float('-inf'), self.__my_true_div(self.x[1], ointerval.x[1])])
            else:
                return Interval([self.__my_true_div(self.x[0], ointerval.x[1]), float('inf')])
        elif ointerval.x[1] == 0:
            if self.x[1] < 0:
                return Interval([self.__my_true_div(self.x[1], ointerval.x[0]), float('inf')])
            elif self.x[0] > 0:
                return Interval([float('-inf'), self.__my_true_div(self.x[0], ointerval.x[0])])
        else: #ointerval.x[0] < 0 < ointerval.x[1]:
            if self.x[1] < 0:
                ra = self.__my_true_div(self.x[1], ointerval.x[1])
                rb = self.__my_true_div(self.x[1], ointerval.x[0])
            elif self.x[0] > 0:
                ra = self.__my_true_div(self.x[0], ointerval.x[0])
                rb = self.__my_true_div(self.x[0], ointerval.x[1])
            return [Interval([float('-inf'), ra]), Interval([rb, float('inf')])]


    def __floordiv__(self, other):
        ointerval = value_to_interval(other)
        if ((not isinstance(self.x[0], int)) and (not isinstance(self.x[0], sympy.Integer))) or ((not isinstance(self.x[1], int)) and
         (not isinstance(self.x[1], sympy.Integer))):
            raise ValueError('floor division is available only for integer bounds.')
        if ((not isinstance(ointerval.x[0], int)) and (not isinstance(ointerval.x[0], sympy.Integer))) or ((not isinstance(ointerval.x[1], int)) and
         (not isinstance(ointerval.x[1], sympy.Integer))):
            raise ValueError('floor division is available only for integer bounds.')
        if ointerval.x[0] == ointerval.x[1] == 0:
            if self.x[0] <= 0 <= self.x[1]:
                return Interval([float('-inf'), float('inf')])
            else:
                return None
        elif ointerval.x[0] > 0 or ointerval.x[1] < 0:
            v = [self.x[0] // ointerval.x[0], self.x[0] // ointerval.x[1], self.x[1] // ointerval.x[0],
                 self.x[1] // ointerval.x[1]]
            b = [min(v), max(v)]
            return Interval(b)
        elif self.x[0] <= 0 <= self.x[1]:
            return Interval([float('-inf'), float('inf')])
        elif ointerval.x[0] == 0:
            if self.x[1] < 0:
                return Interval([float('-inf'), self.__my_floor_div(self.x[1], ointerval.x[1])])
            else:
                return Interval([self.__my_floor_div(self.x[0], ointerval.x[1]), float('inf')])
        elif ointerval.x[1] == 0:
            if self.x[1] < 0:
                return Interval([self.__my_floor_div(self.x[1], ointerval.x[0]), float('inf')])
            elif self.x[0] > 0:
                return Interval([float('-inf'), self.__my_floor_div(self.x[0], ointerval.x[0])])
        else: #ointerval.x[0] < 0 < ointerval.x[1]:
            if self.x[1] < 0:
                ra = self.__my_floor_div(self.x[1], ointerval.x[1])
                rb = self.__my_floor_div(self.x[1], ointerval.x[0])
            elif self.x[0] > 0:
                ra = self.__my_floor_div(self.x[0], ointerval.x[0])
                rb = self.__my_floor_div(self.x[0], ointerval.x[1])
            return [Interval([float('-inf'), ra]), Interval([rb, float('inf')])]

    def __rmul__(self, other):
        return self.__mul__(other)

    def __rtruediv__(self, other):
        ointerval = value_to_interval(other)
        return ointerval.__truediv__(self)

    def __rfloordiv__(self, other):
        ointerval = value_to_interval(other)
        return ointerval.__floordiv__(self)

    def is_dot_interval(self):
        """
        Checks whether the given interval is a dot interval, i.e. left end = right end

        Returns:
        --------
        True if x is a dot interval, False otherwise
        """
        return self.x[0] == self.x[1]

class IntervalMP(Interval):
    x = [None, None]
    def __init__(self, x):
        """
        Constructor for Interval

        Parameters
        ----------
        x - list with size 2 [a, b]. Where a < b. a - lower bound, b - upper bound. a, b: Decimal or Interval
        """
        if isinstance(x, Interval): #Добавить в ТЗ!!!
            self.x[0] = dec.Decimal(x.x[0])
            self.x[1] = dec.Decimal(x.x[1])
            return
        if len(x) != 2:
            raise ValueError("Length of the list in parameters must be 2.")
        if x[0] > x[1]:
            return None
        if type(x[0]) != dec.Decimal or type(x[1]) != dec.Decimal:
             raise TypeError("Interval constructor's arguments must be instances of Decimal")

        self.x = x.copy()


    def mid(self):
        """
        Returns the best approximation of the interval central point
        """
        prev = set_rounding_mode_default()
        mid =  dec.Decimal('0.5') * (self.x[0] + self.x[1])
        set_rounding_mode(prev)
        return mid


    def width(self, rounding_mode = dec.getcontext().rounding):
        """
        Returns the closest interval's width with given approximation mode.

        Parameters
        ----------
        rounding_mode - rounding mode for approximation from decimal library.
                        Current rounding by default.
        """
        prev = set_rounding_mode(rounding_mode)
        width = self.x[1] - self.x[0]
        set_rounding_mode(prev)
        return width

    def radius(self, rounding_mode = dec.getcontext().rounding):
        """
        Returns the closest interval's radius with given approximation mode.

        Parameters
        ----------
        rounding_mode - rounding mode for approximation from decimal library.
                        Current rounding by default.
        """
        prev = set_rounding_mode(rounding_mode)
        radius = dec.Decimal('0.5') * (self.x[0] - self.x[1])
        set_rounding_mode(prev)
        return radius

    def scale(self, factor, rounding_mode = dec.getcontext().rounding):
        """
        Changes width of the interval in given number of times with respect to the middle.
        If rounding_mode = dec.ROUND_ROUND_HALF_EVEN, calculetes the closest approximation of the new bounds.
        If rounding_mode = dec.ROUND_CEILING, calculetes the approximation of the new bounds, when width is it's maximum.
        If rounding_mode = dec.ROUND_FLOOR, calculetes the approximation of the new bounds, when width is it's minimum.


        Parameters
        ----------
        factor - number of times in which interval's width is changing.
        rounding_mode - rounding mode for approximation from decimal library.
                        Current rounding by default.
        """
        prev = dec.getcontext().rounding
        m = self.mid()
        r = self.radius(rounding_mode)
        new_r = r * factor
        self.x[1] = m + new_r
        if rounding_mode == dec.ROUND_CEILING:
            set_rounding_mode_floor()
        elif rounding_mode == dec.ROUND_FLOOR:
            set_rounding_mode_ceil()
        self.x[0] = m - new_r
        set_rounding_mode(prev)

    def intersec(self, other):
        """
        Returns intersection of two intervals

        Parameters
        ----------
        other - interval or value
        """
        ointerval = value_to_interval_mp(other)
        if self.x[1] < ointerval.x[0] or ointerval.x[1] < self.x[0]:
            return None
        return IntervalMP([max(self.x[0], ointerval.x[0]), min(self.x[1], ointerval.x[1])])

    def __neg__(self):
        ninterval = IntervalMP(self.x)
        ninterval.x[0] = - self.x[1]
        ninterval.x[1] = - self.x[0]
        return ninterval

    def __add__(self, other):
        prev = dec.getcontext().rounding
        ointerval = value_to_interval_mp(other)
        ninterval = IntervalMP(self.x)
        ninterval.x[1] = self.x[1] + ointerval.x[1]
        if prev == dec.ROUND_CEILING:
            set_rounding_mode_floor()
        elif prev == dec.ROUND_FLOOR:
            set_rounding_mode_ceil()
        ninterval.x[0] = self.x[0] + ointerval.x[0]
        set_rounding_mode(prev)
        return ninterval

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        ointerval = value_to_interval_mp(other)
        return self.__add__(-ointerval)

    def __rsub__(self, other):
        ointerval = value_to_interval_mp(other)
        return (-self).__add__(ointerval)


    def __pow__(self, factor):
        if ((not isinstance(factor, int)) and (not isinstance(factor, sympy.Integer))) and ((not isinstance(factor, float)) and
         (not isinstance(factor, sympy.Float))) and (not isinstance(factor, dec.Decimal)):
            raise TypeError("Factor for exponentiation should be int, float or decimal")
        is_nec_fact = False
        if factor < 0:
          factor *= -1
          is_nec_fact = True
        ninterval = IntervalMP(self.x)
        prev = dec.getcontext().rounding
        set_rounding_mode_ceil()
        u_ceil = self.x[0] ** factor
        v_ceil = self.x[1] ** factor
        set_rounding_mode_default()
        u_default = self.x[0] ** factor
        v_default = self.x[1] ** factor
        set_rounding_mode_floor()
        u_floor = self.x[0] ** factor
        v_floor = self.x[1] ** factor
        if factor == 0:
            ninterval.x[0] = dec.Decimal('1')
            ninterval.x[1] = dec.Decimal('1')
        elif factor % 2 == 0:
            if prev == dec.ROUND_HALF_EVEN:
                ninterval.x[1] = max(u_default, v_default)
                if self.x[0] <= dec.Decimal('0') and self.x[1] >= dec.Decimal('0'):
                    ninterval.x[0] = dec.Decimal('0')
                else:
                    ninterval.x[0] = min(u_default, v_default)
            elif prev == dec.ROUND_CEILING:
                ninterval.x[1] = max(u_ceil, v_ceil)
                if self.x[0] <= dec.Decimal('0') and self.x[1] >= dec.Decimal('0'):
                    ninterval.x[0] = dec.Decimal('0')
                else:
                    ninterval.x[0] = min(u_floor, v_floor)
            else:
                ninterval.x[1] = max(u_floor, v_floor)
                if self.x[0] <= dec.Decimal('0') and self.x[1] >= dec.Decimal('0'):
                    ninterval.x[0] = dec.Decimal('0')
                else:
                    ninterval.x[0] = min(u_ceil, v_ceil)
        else:
            if prev == dec.ROUND_HALF_EVEN:
                ninterval.x[0] = min(u_default, v_default)
                ninterval.x[1] = max(u_default, v_default)
            elif prev == dec.ROUND_CEILING:
                ninterval.x[0] = min(u_floor, v_floor)
                ninterval.x[1] = max(u_ceil, v_ceil)
            else:
                ninterval.x[0] = min(u_ceil, v_ceil)
                ninterval.x[1] = max(u_floor, v_floor)
        if is_nec_fact:
            ninterval = 1/ninterval
        set_rounding_mode(prev)
        return ninterval

    def __rpow__(self, value):
        if ((not isinstance(value, int)) and (not isinstance(value, sympy.Integer))) and ((not isinstance(value, float)) and
         (not isinstance(value, sympy.Float))) and (not isinstance(value, dec.Decimal)):
            raise TypeError("Value for exponentiation should be int, float or decimal")
        if value < 0:
            raise ValueError("Value should be positive or 0")
        if value == 0:
            if self.x[0] < dec.Decimal('0'):
                raise ValueError("0 can't be in negative power")
            elif self.x[0] == dec.Decimal('0') and self.x[1] == dec.Decimal('0'):
                return IntervalMP([dec.Decimal('1'), dec.Decimal('1')])
            elif self.x[0] == dec.Decimal('0') and self.x[1] > dec.Decimal('0'):
                return [IntervalMP([dec.Decimal('1'), dec.Decimal('1')]), IntervalMP([dec.Decimal('0'), dec.Decimal('0')])]
            else:
                prev = dec.getcontext().rounding
                if prev == dec.ROUND_HALF_EVEN:
                    u = value**self.x[0]
                    v = value**self.x[1]
                    return IntervalMP([min(u, v), max(u, v)])

                set_rounding_mode_ceil()
                u_ceil = value**self.x[0]
                v_ceil = value**self.x[1]
                set_rounding_mode_floor()
                u_floor = value**self.x[0]
                v_floor = value**self.x[1]
                set_rounding_mode(prev)
                if prev == dec.ROUND_CEILING:
                    return IntervalMP([min(u_floor, v_floor), max(u_ceil, v_ceil)])
                else:
                    return IntervalMP([min(u_ceil, v_ceil), max(u_floor, v_floor)])

    def __mul__(self, other):
        ointerval = value_to_interval_mp(other)
        prev = dec.getcontext().rounding
        if prev == dec.ROUND_HALF_EVEN:
            v = [self.x[0] * ointerval.x[0], self.x[0] * ointerval.x[1], self.x[1] * ointerval.x[0], self.x[1] * ointerval.x[1]]
            b = [min(v), max(v)]
        set_rounding_mode_ceil()
        v_ceil = [self.x[0] * ointerval.x[0], self.x[0] * ointerval.x[1], self.x[1] * ointerval.x[0], self.x[1] * ointerval.x[1]]
        set_rounding_mode_floor()
        v_floor = [self.x[0] * ointerval.x[0], self.x[0] * ointerval.x[1], self.x[1] * ointerval.x[0], self.x[1] * ointerval.x[1]]
        set_rounding_mode(prev)
        if prev == dec.ROUND_CEILING:
            b = [min(v_floor), max(v_ceil)]
        else:
            b = [min(v_ceil), max(v_floor)]
        return IntervalMP(b)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        ointerval = value_to_interval_mp(other)
        prev = dec.getcontext().rounding
        if ointerval.x[0] == ointerval.x[1] == dec.Decimal('0'):
            if self.x[0] <= dec.Decimal('0') <= self.x[1]:
                return IntervalMP([dec.Decimal('-Infinity'), dec.Decimal('Infinity')])
            else:
                return None
        elif ointerval.x[0] > dec.Decimal('0') or ointerval.x[1] < dec.Decimal('0'):
            u = 1 / ointerval.x[0]
            if prev == dec.ROUND_CEILING:
                set_rounding_mode_floor()
            elif prev == dec.ROUND_FLOOR:
                set_rounding_mode_ceil()
            v = 1 / ointerval.x[1]
            set_rounding_mode(prev)
            return self.__mul__(IntervalMP([v, u]))
        elif self.x[0] <= dec.Decimal('0') <= self.x[1]:
            return Interval([dec.Decimal('-Infinity'), dec.Decimal('Infinity')])
        elif ointerval.x[0] == dec.Decimal('0'):
            if self.x[1] < dec.Decimal('0'):
                u = self.__my_true_div(self.x[1], ointerval.x[1])
                return Interval([dec.Decimal('-Infinity'), u])
            else:
                if prev == dec.ROUND_CEILING:
                    set_rounding_mode_floor()
                elif prev == dec.ROUND_FLOOR:
                    set_rounding_mode_ceil()
                v = self.__my_true_div(self.x[0], ointerval.x[1])
                set_rounding_mode(prev)
                return Interval([v, dec.Decimal('Infinity')])
        elif ointerval.x[1] == dec.Decimal('0'):
            if self.x[1] < dec.Decimal('0'):
                if prev == dec.ROUND_CEILING:
                    set_rounding_mode_floor()
                elif prev == dec.ROUND_FLOOR:
                    set_rounding_mode_ceil()
                v = self.__my_true_div(self.x[1], ointerval.x[0])
                set_rounding_mode(prev)
                return Interval([v, dec.Decimal('Infinity')])
            elif self.x[0] > dec.Decimal('0'):
                u = self.__my_true_div(self.x[0], ointerval.x[0])
                return Interval([dec.Decimal('-Infinity'), u])
        else:
            if self.x[1] < dec.Decimal('0'):
                ra = self.__my_true_div(self.x[1], ointerval.x[1])
                if prev == dec.ROUND_CEILING:
                    set_rounding_mode_floor()
                elif prev == dec.ROUND_FLOOR:
                    set_rounding_mode_ceil()
                rb = self.__my_true_div(self.x[1], ointerval.x[0])
            elif self.x[0] > dec.Decimal('0'):
                ra = self.__my_true_div(self.x[0], ointerval.x[0])
                if prev == dec.ROUND_CEILING:
                    set_rounding_mode_floor()
                elif prev == dec.ROUND_FLOOR:
                    set_rounding_mode_ceil()
                rb = self.__my_true_div(self.x[0], ointerval.x[1])
            set_rounding_mode(prev)
            return [IntervalMP([dec.Decimal('-Infinity'), ra]), Interval([rb, dec.Decimal('Infinity')])]

    def __rtruediv__(self, other):
        ointerval = value_to_interval(other)
        return ointerval.__truediv__(self)

    def is_dot_interval(self):
        """
        Checks whether the given interval is a dot interval, i.e. left end = right end

        Returns:
        --------
        True if x is a dot interval, False otherwise
        """
        res = dec.compare(self.a, self.b)
        return res == dec.Decimal('0')



def set_precision(prec):
    """
    Sets the precision as the number of significant decimal digits in mantissa for decimal values.

    Parameters
    ----------
    prec : integer
          The number of decimal places (should be between 1 and some reasonable value)
    """

    prev = dec.getcontext().prec
    dec.getcontext().prec = prec
    return prev

def set_rounding_mode(rounding_mode):
    """
    Changes current rouding mode for decimal values.

    Parameters
    ----------
    rounding_mode - rounding mode for approximation from decimal library.
    """
    prev = dec.getcontext().rounding
    dec.getcontext().rounding = rounding_mode
    return prev


def set_rounding_mode_default():
    """
    Changes current rouding mode to ROUND_HALF_EVEN for decimal values.
    """
    return set_rounding_mode(dec.ROUND_HALF_EVEN)


def set_rounding_mode_ceil():
    """
    Changes current rouding mode to ROUND_CEILING for decimal values.
    """
    return set_rounding_mode(dec.ROUND_CEILING)


def set_rounding_mode_floor():
    """
    Changes current rouding mode to ROUND_FLOOR for decimal values.
    """
    return set_rounding_mode(dec.ROUND_FLOOR)


def value_to_interval(expr):
    """
    Returns Interval, which equals the parameter.

    Parameters
    ----------
    expr - int, float, sympy.Integer, sympy.Float, dec.Decimal or Interval.

    """
    if isinstance(expr, int):
        etmp = Interval([expr, expr])
    elif isinstance(expr, float):
        etmp = Interval([expr, expr])
    elif isinstance(expr, sympy.Integer):
        etmp = Interval([int(expr), int(expr)])
    elif isinstance(expr, (sympy.Float, dec.Decimal)):
        etmp = Interval([float(expr), float(expr)])
    elif  isinstance(expr, Interval):
        etmp = expr
    else:
        raise TypeError("The parameter should be a number or an Interval.")
    return etmp

def value_to_interval_mp(expr):
    """
    Returns IntervalMP, which equals the parameter.

    Parameters
    ----------
    expr - int, float, sympy.Integer, sympy.Float or Interval.
    """
    if isinstance(expr, int):
        etmp = IntervalMP([dec.Decimal(expr), dec.Decimal(expr)])
    elif isinstance(expr, (float, dec.Decimal)):
        etmp = IntervalMP([dec.Decimal(expr), dec.Decimal(expr)])
    elif isinstance(expr, sympy.Integer):
        etmp = IntervalMP([dec.Decimal(int(expr)), dec.Decimal(int(expr))])
    elif isinstance(expr, sympy.Float):
        etmp = IntervalMP([dec.Decimal(float(expr)), dec.Decimal(float(expr))])
    elif isinstance(expr, IntervalMP):
        etmp = expr
    elif isinstance(expr, Interval) and not isinstance(expr, IntervalMP):
        etmp = IntervalMP(expr)
    else:
        raise TypeError("The parameter should be a number or an Interval.")
    return etmp


def _pi():
    """
    Compute decimal Pi to the current precision.
    """
    dec.getcontext().prec += 2
    three = dec.Decimal(3)
    lasts = 0
    t = three
    s = 3
    n = 1
    na = 0
    d = 0
    da = 24
    while s != lasts:
        lasts = s
        n = n + na
        na += 8
        d = d + da
        da += 32
        t = (t * n) / d
        s += t
    dec.getcontext().prec -= 2
    return +s


def _sin_dec(x):
    """
    Compute decimal sin(x) to the current precision.

    Parameters
    ----------
    x - Decimal
    """
    pi_dec = _pi()
    x = x % (2 * pi_dec)
    dec.getcontext().prec += 2
    i = 1
    lasts = 0
    s = x
    fact = 1
    num = x
    sign = 1
    while s != lasts:
        lasts = s
        i += 2
        fact *= i * (i-1)
        num *= x * x
        sign *= -1
        s += num / fact * sign
    dec.getcontext().prec -= 2
    return +s

def _cos_dec(x):
    """
    Compute decimal cos(x) to the current precision.

    Parameters
    ----------
    x - Decimal
    """
    pi_dec = _pi()
    x = x % (2 * pi_dec)
    dec.getcontext().prec += 2
    i = 0
    lasts = 0
    s = 1
    fact = 1
    num = 1
    sign = 1
    while s != lasts:
        lasts = s
        i += 2
        fact *= i * (i-1)
        num *= x * x
        sign *= -1
        s += num / fact * sign
    dec.getcontext().prec -= 2
    return +s

def _sin(x):
    """
    Compute sin(x).

    Parameters
    ----------
    x - int, float, Decimal, Interval
    """
    if isinstance(x, (int, sympy.Integer, np.integer)):
        return math.sin(x)
    elif isinstance(x, (float, sympy.Float, np.float32, np.float64)):
        return math.sin(x)
    elif isinstance(x, dec.Decimal):
        return _sin_dec(x)
    else:
        if isinstance(x, IntervalMP):
            prev = set_rounding_mode_default()
            pi_dec = _pi()
            pi2 =  dec.Decimal('2') * pi_dec
            pi05 = pi_dec * dec.Decimal('0.5')
            y_def = [_sin_dec(x[0]), _sin_dec(x[1])]
            set_rounding_mode_ceil()
            y_ceil = [_sin_dec(x[0]), _sin_dec(x[1])]
            set_rounding_mode_floor()
            y_floor = [_sin_dec(x[0]), _sin_dec(x[1])]
            set_rounding_mode(prev)
            if math.ceil(float((x[0] - pi05)/pi2)) <= math.floor((x[1] - pi05)/pi2):
                b = dec.Decimal('1')
            elif prev == dec.ROUND_HALF_EVEN:
                b = max(y_def)
            elif prev == dec.ROUND_FLOOR:
                b = max(y_floor)
            else:
                b = max(y_ceil)
            if math.ceil(float((x[0] + pi05)/pi2)) <= math.floor(float((x[1] + pi05)/pi2)):
                a = dec.Decimal('-1')
            elif prev == dec.ROUND_HALF_EVEN:
                a = min(y_def)
            elif prev == dec.ROUND_FLOOR:
                a = min(y_ceil)
            else:
                a = min(y_floor)
            return IntervalMP([a, b])
        y = [math.sin(x[0]), math.sin(x[1])]
        pi2 = 2 * math.pi
        pi05 = math.pi / 2
        if math.ceil((x[0] - pi05)/pi2) <= math.floor((x[1] - pi05)/pi2):
            b = 1
        else:
            b = max(y)
        if math.ceil((x[0] + pi05)/pi2) <= math.floor((x[1] + pi05)/pi2):

            a = -1
        else:
            a = min(y)
        return Interval([a,b])

def sin(x):
    """
    Compute sin(x).

    Parameters
    ----------
    x - int, float, Decimal, Interval or list or arrays of ints, floats, Decimals or Intervals.
    """
    if isinstance(x, (int, sympy.Integer, np.integer)):
        return math.sin(x)
    elif isinstance(x, (float, sympy.Float, np.float32, np.float64)):
        return math.sin(x)
    elif isinstance(x, dec.Decimal):
        return _sin_dec(x)
    elif isinstance(x, Interval):
        return _sin(x)
    x = np.array(x)
    buf = x.ravel()
    outvector = np.empty_like(buf)
    for i in range(buf.shape[0]):
        outvector[i] = _sin(buf[i])
    return outvector.reshape(x.shape)

def _cos(x):
    """
    Compute cos(x).

    Parameters
    ----------
    x - int, float, Decimal, Interval
    """
    if isinstance(x, (int, sympy.Integer, np.integer)):
        return math.cos(x)
    elif isinstance(x, (float, sympy.Float, np.float32, np.float64)):
        return math.cos(x)
    elif isinstance(x, dec.Decimal):
        return _cos_dec(x)
    else:
        if isinstance(x, IntervalMP):
            prev = set_rounding_mode_default()
            pi_dec = _pi()
            pi2 =  dec.Decimal('2') * pi_dec
            y_def = [_cos_dec(x[0]), _cos_dec(x[1])]
            set_rounding_mode_ceil()
            y_ceil = [_cos_dec(x[0]), _cos_dec(x[1])]
            set_rounding_mode_floor()
            y_floor = [_cos_dec(x[0]), _cos_dec(x[1])]
            set_rounding_mode(prev)
            if math.ceil(float(x[0]/pi2)) <= math.floor(float(x[1]/pi2)):
                b = dec.Decimal('1')
            elif prev == dec.ROUND_HALF_EVEN:
                b = max(y_def)
            elif prev == dec.ROUND_FLOOR:
                b = max(y_floor)
            else:
                b = max(y_ceil)

            if math.ceil(float((x[0] - pi_dec)/pi2)) <=  math.floor(float((x[1] - pi_dec)/pi2)):

                a = dec.Decimal('-1')
            elif prev == dec.ROUND_HALF_EVEN:
                a = min(y_def)
            elif prev == dec.ROUND_FLOOR:
                a = min(y_ceil)
            else:
                a = min(y_floor)
            return IntervalMP([a, b])
        y = [math.cos(x[0]), math.cos(x[1])]
        pi2 = 2 * math.pi
        if math.ceil(x[0]/pi2) <= math.floor(x[1]/pi2):
            b = 1
        else:
            b = max(y)
        if math.ceil((x[0] - math.pi)/pi2) <= math.floor((x[1] - math.pi)/pi2):
            a = -1
        else:
            a = min(y)
        return Interval([a,b])

def cos(x):
    """
    Compute cos(x).

    Parameters
    ----------
    x - int, float, Decimal, Interval or list or arrays of ints, floats, Decimals or Intervals.
    """
    if isinstance(x, (int, sympy.Integer, np.integer)):
        return math.cos(x)
    elif isinstance(x, (float, sympy.Float, np.float32, np.float64)):
        return math.cos(x)
    elif isinstance(x, dec.Decimal):
        return _cos_dec(x)
    elif isinstance(x, Interval):
        return _cos(x)
    x = np.array(x)
    buf = x.ravel()
    outvector = np.empty_like(buf)
    for i in range(buf.shape[0]):
        outvector[i] = _cos(buf[i])
    return outvector.reshape(x.shape)


def _exp(x):
    """
    Compute e^x.

    Parameters
    ----------
    x - int, float, Decimal, Interval
    """
    if isinstance(x, (int, np.integer)):
        return math.exp(x)
    elif isinstance(x, (float, np.float32, np.float64)):
        return math.exp(x)
    elif isinstance(x, dec.Decimal):
        return x.exp()
    else:
        if isinstance(x, IntervalMP):
            u = x[1].exp()
            v = x[0].exp()
            prev = dec.getcontext().rounding
            if prev == dec.ROUND_FLOOR:
                set_rounding_mode_ceil()
                v = x[0].exp()
            elif prev == dec.ROUND_CEILING:
                set_rounding_mode_floor()
                v = x[0].exp()
            set_rounding_mode(prev)
            return IntervalMP([v, u])
        return Interval([math.exp(x[0]), math.exp(x[1])])

def exp(x):
    """
    Compute e^x.

    Parameters
    ----------
    x - int, float, Decimal, Interval or list/array/matrix of ints, floats, Decimals or Intervals.
    """
    if isinstance(x, (int, np.integer)):
        return math.exp(x)
    elif isinstance(x, (float, np.float32, np.float64)):
        return math.exp(x)
    elif isinstance(x, dec.Decimal):
        return x.exp()
    elif isinstance(x, Interval):
        return _exp(x)
    x = np.array(x)
    buf = x.ravel()
    outvector = np.empty_like(buf)
    for i in range(buf.shape[0]):
        outvector[i] = _exp(buf[i])
    return outvector.reshape(x.shape)

def _abs(x):
    """
    Returns absolute value of x, where x is number or Interval

    """
    if isinstance(x, (int, np.integer, float, np.float32, np.float64)):
        if x < 0:
            return -x
        return x
    elif isinstance(x, dec.Decimal):
        if x < dec.Decimal('0'):
            return -x
        return x
    elif isinstance(x, IntervalMP):
        if x[1] <= dec.Decimal('0'):
            return IntervalMP([-x[1], -x[0]])
        elif x[0] <= dec.Decimal('0') and x[1] > dec.Decimal('0'):
            if -x[0] > x[1]:
              return IntervalMP([dec.Decimal('0'), -x[0]])
            else:
              return IntervalMP([dec.Decimal('0'), x[1]])
        else:
            return IntervalMP([x[0], x[1]])
    elif isinstance(x, Interval):
        if x[1] <= 0:
            return Interval([-x[1], -x[0]])
        elif x[0] <= 0 and x[1] > 0:
            if -x[0] > x[1]:
              return Interval([0, -x[0]])
            else:
              return Interval([0, x[1]])
        else:
            return Interval([x[0], x[1]])

def abs(x):
    """
    Returns absolute value of x, where x is number, Interval or list/array/ matrix of numbers or Intervals.
    """
    if isinstance(x, (int, np.integer, float, np.float32, np.float64, dec.Decimal)):
        return _abs(x)
    elif isinstance(x, Interval):
        return _abs(x)
    x = np.array(x)
    buf = x.ravel()
    outvector = np.empty_like(buf)
    for i in range(buf.shape[0]):
        outvector[i] = _abs(buf[i])
    return outvector.reshape(x.shape)

def _log_dec(x, base):
    """
    Compute decimal logarithm of x to base to the current precision, where at least one of the parameters is Decimal.
    and other is a number - float or int.
    """
    if not isinstance(x, dec.Decimal) and isinstance(base, dec.Decimal):
        return dec.Decimal(math.log(x))/base.ln()
    elif isinstance(x, dec.Decimal) and not isinstance(base, dec.Decimal):
        return x.ln()/dec.Decimal(math.log(base))
    elif isinstance(x, dec.Decimal) and isinstance(base, dec.Decimal):
        return x.ln()/base.ln()

def _log(x, base):
    """
    Compute logarithm of x to base, where x could be an Interval or a number.

    Parameters
    ----------
    x - int, float, Decimal or Interval.
    base - int, float or Decimal.
    """
    if not isinstance(base, (int, np.integer, sympy.Integer, float, np.float32, np.float64, sympy.Float, dec.Decimal)):
        raise TypeError("Base should be int, float or decimal")
    if isinstance(x, (int, np.integer, sympy.Integer, float, np.float32, np.float64, sympy.Float)) and\
    isinstance(base, (int, np.integer, sympy.Integer, float, np.float32, np.float64, sympy.Float)):
        return math.log(x, base)
    elif isinstance(x, dec.Decimal):
        return _log_dec(x, base)
    elif isinstance(x, (int, np.integer, sympy.Integer, float, np.float32, np.float64, sympy.Float)) and isinstance(base, dec.Decimal):
        return _log_dec(x, base)
    elif isinstance(x, IntervalMP):
        prev = set_rounding_mode_default()
        y_def = [_log_dec(x[0], base), _log_dec(x[1], base)]
        set_rounding_mode_ceil()
        y_ceil = [_log_dec(x[0], base), _log_dec(x[1], base)]
        set_rounding_mode_floor()
        y_floor = [_log_dec(x[0], base), _log_dec(x[1], base)]
        set_rounding_mode(prev)
        if dec.Decimal(base) > dec.Decimal('1'):
            if prev == dec.ROUND_CEILING:
                return IntervalMP([y_floor[0], y_ceil[1]])
            elif prev == dec.ROUND_FLOOR:
                return IntervalMP(y_ceil[0], y_floor[1])
            return IntervalMP([y_def])
        else:
            if prev == dec.ROUND_CEILING:
                return IntervalMP([y_floor[1], y_ceil[0]])
            elif prev == dec.ROUND_FLOOR:
                return IntervalMP([y_ceil[1], y_floor[0]])
            return IntervalMP([y_def[1], y_def[0]])
    elif isinstance(x, Interval):
        if base > 1:
            return Interval([math.log(x[0], base), math.log(x[1], base)])
        else:
            return Interval([math.log(x[1], base), math.log(x[0], base)])


def log(x, base):
    """
    Compute logarithm of x to base, where x could be an Interval, a number or list/array/matrix.

    Parameters
    ----------
    x - int, float, Decimal, Interval or list/array/matrix of ints, floats, Decimals or Intervals.
    base - int, float or Decimal.
    """
    if not isinstance(base, (int, np.integer,sympy.Integer, float, np.float32, np.float64, sympy.Float, dec.Decimal)):
        raise TypeError("Base should be int float or decimal")
    if isinstance(x, (int, np.integer, float, sympy.Integer, np.float32, np.float64, sympy.Float, dec.Decimal, Interval)):
        return _log(x, base)
    x = np.array(x)
    buf = x.ravel()
    outvector = np.empty_like(buf)
    for i in range(buf.shape[0]):
        outvector[i] = _log(buf[i], base)
    return outvector.reshape(x.shape)


def diff(ival, fun, x):
    """
    Returns an Interval - the first derivative fun from x, in the interval ival.
    The return type is always Interval without machine precision.

    Parameters
    ----------
    ival - int, float, Decimal, Interval.
    fun - function with parameter x.
    x - sympy.Symbol, which is the parameter of fun.
    """
    ival = value_to_interval(ival)
    if not isinstance(x, sympy.Symbol):
        raise TypeError("Function parameters must be ival: Interval, fun: function, x: sympy.Symbol")
    diff_1 = sympy.diff(fun, x)
    diff_2 = sympy.diff(diff_1, x)
    if diff_2 == 0:
        return Interval([float(diff_1),float(diff_1)])
    roots = sympy.solveset(diff_2, domain = sympy.S.Reals)
    up = max(diff_1.subs(x, ival[0]), diff_1.subs(x, ival[1]))
    down = min(diff_1.subs(x, ival[0]), diff_1.subs(x, ival[1]))
    if len(roots) == 0:
        return Interval([down, up])
    for i in roots:
        buf = value_to_interval(float(i))
        if buf.is_in(ival):
            up = max(up, diff_1.subs(x, i))
            down = min(down, diff_1.subs(x, i))
    return Interval([down, up])


def mid(x):
    """
    Returns the middle of x.
    If x is a number returns x.
    If x is an Interval returns the middle of the interval.
    If x is a list/array/matrix of ints, floats, Decimals or Intervals returns np.array with middles for each value.
    """
    if isinstance(x, (int, np.integer,sympy.Integer, float, np.float32, np.float64, sympy.Float, dec.Decimal)):
        return x
    if isinstance(x, Interval):
        return x.mid()
    x = np.array(x)
    buf = x.ravel()
    outvector = np.empty_like(buf)
    for i in range(buf.shape[0]):
        outvector[i] = buf[i].mid()
    return outvector.reshape(x.shape)

def width(x):
    """
    Returns the width of x.
    If x is a number returns 0.
    If x is an Interval returns the width of the interval.
    If x is a list/array/matrix of ints, floats, Decimals or Intervals returns np.array with widths for each value.
    """
    if isinstance(x, (int, np.integer,sympy.Integer, float, np.float32, np.float64, sympy.Float, dec.Decimal)):
        return 0
    if isinstance(x, Interval):
        return x.width()
    x = np.array(x)
    buf = x.ravel()
    outvector = np.empty_like(buf)
    for i in range(buf.shape[0]):
        outvector[i] = buf[i].width()
    return outvector.reshape(x.shape)

def radius(x):
    """
    Returns the radius of x.
    If x is a number returns 0.
    If x is an Interval returns the radius of the interval.
    If x is a list/array/matrix of ints, floats, Decimals or Intervals returns np.array with radiuses for each value.
    """
    if isinstance(x, (int, np.integer,sympy.Integer, float, np.float32, np.float64, sympy.Float, dec.Decimal)):
        return 0
    if isinstance(x, Interval):
        return x.radius()
    x = np.array(x)
    buf = x.ravel()
    outvector = np.empty_like(buf)
    for i in range(buf.shape[0]):
        outvector[i] = buf[i].radius()
    return outvector.reshape(x.shape)