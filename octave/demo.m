sys = tds_create({1.0, -1.0}, [0.0, 1.0]);

l = tds_roots(sys, 10);

plot(real(l), imag(l));