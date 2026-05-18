package main

import (
	"fmt"
	"math"
)

type Complex struct {
	Re, Im float64
}

func NewComplex(re, im float64) Complex {
	return Complex{Re: re, Im: im}
}

func (c Complex) Add(other Complex) Complex {
	return Complex{Re: c.Re + other.Re, Im: c.Im + other.Im}
}

func (c Complex) Mul(other Complex) Complex {
	return Complex{
		Re: c.Re*other.Re - c.Im*other.Im,
		Im: c.Re*other.Im + c.Im*other.Re,
	}
}

func (c Complex) MulScalar(r float64) Complex {
	return Complex{Re: c.Re * r, Im: c.Im * r}
}

func (c Complex) String() string {
	if math.Abs(c.Im) < 1e-12 {
		return fmt.Sprintf("%.6f", c.Re)
	}
	if math.Abs(c.Re) < 1e-12 {
		return fmt.Sprintf("%.6fi", c.Im)
	}
	sign := '+'
	if c.Im < 0 {
		sign = '-'
	}
	return fmt.Sprintf("%.6f %c %.6fi", c.Re, sign, math.Abs(c.Im))
}

func cbrtComplex(z Complex) Complex {
	r := math.Pow(math.Hypot(z.Re, z.Im), 1.0/3.0)
	theta := math.Atan2(z.Im, z.Re) / 3.0
	return Complex{
		Re: r * math.Cos(theta),
		Im: r * math.Sin(theta),
	}
}

func SolveCubic(a, b, c, d float64) [3]Complex {
	if math.Abs(a) < 1e-14 {
		panic("a is zero, not a cubic equation")
	}

	p := b / a
	q := c / a
	r := d / a

	P := q - p*p/3.0
	Q := (2.0*p*p*p)/27.0 - (p*q)/3.0 + r

	Delta := (Q*Q)/4.0 + (P*P*P)/27.0

	if math.Abs(Delta) < 1e-12 {
		Delta = 0
	}

	var t [3]Complex

	if Delta > 0 {
		sqrtDelta := math.Sqrt(Delta)
		uReal := math.Cbrt(-Q/2.0 + sqrtDelta)
		vReal := math.Cbrt(-Q/2.0 - sqrtDelta)

		u := NewComplex(uReal, 0)
		v := NewComplex(vReal, 0)

		t[0] = u.Add(v)

		omega := NewComplex(-0.5, math.Sqrt(3)/2.0)
		omega2 := NewComplex(-0.5, -math.Sqrt(3)/2.0)

		t[1] = u.Mul(omega).Add(v.Mul(omega2))
		t[2] = u.Mul(omega2).Add(v.Mul(omega))
	} else if math.Abs(Delta) < 1e-12 {
		uReal := math.Cbrt(-Q / 2.0)
		t[0] = NewComplex(2.0*uReal, 0)
		t[1] = NewComplex(-uReal, 0)
		t[2] = NewComplex(-uReal, 0)
	} else {
		rho := math.Sqrt(-P*P*P/27.0) * 2.0
		phi := math.Acos((-Q / 2.0) / math.Sqrt(-P*P*P/27.0))

		t[0] = NewComplex(rho*math.Cos(phi/3.0), 0)
		t[1] = NewComplex(rho*math.Cos((phi+2.0*math.Pi)/3.0), 0)
		t[2] = NewComplex(rho*math.Cos((phi+4.0*math.Pi)/3.0), 0)
	}

	for i := 0; i < 3; i++ {
		t[i] = t[i].Add(NewComplex(-p/3.0, 0))
	}

	return t
}

func prettyPrint(roots [3]Complex) {
	for i, root := range roots {
		fmt.Printf("x%d = %s\n", i+1, root.String())
	}
}

func main() {
	fmt.Println("Решение кубического уравнения a*x^3 + b*x^2 + c*x + d = 0")
	fmt.Println("Формула Кардано (с комплексными корнями)\n")

	fmt.Println("Пример 1: 1*x^3 - 6*x^2 + 11*x - 6 = 0")
	roots1 := SolveCubic(1, -6, 11, -6)
	prettyPrint(roots1)
	fmt.Println()

	fmt.Println("Пример 2: 1*x^3 + 0*x^2 + 1*x + 0 = 0")
	roots2 := SolveCubic(1, 0, 1, 0)
	prettyPrint(roots2)
	fmt.Println()

	fmt.Println("Пример 3: 1*x^3 - 3*x^2 + 3*x - 1 = 0")
	roots3 := SolveCubic(1, -3, 3, -1)
	prettyPrint(roots3)
	fmt.Println()

	fmt.Println("Пример 4: 1*x^3 + 0*x^2 - 2*x - 4 = 0")
	roots4 := SolveCubic(1, 0, -2, -4)
	prettyPrint(roots4)
	fmt.Println()

	fmt.Println("Пример 5: 2*x^3 - 4*x^2 - 22*x + 24 = 0")
	roots5 := SolveCubic(2, -4, -22, 24)
	prettyPrint(roots5)
}
