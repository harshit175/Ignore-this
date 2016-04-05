syms  A a11 a12 a13 a14 a21 a22 a23 a24 a31 a32 a33 a34 a41 a42 a43 a44 B b11 b22 b33 b44 C c11 c22 c33 c44 D d11 d22 d33 d44 real

%A=[a11 a12 a13 a14; a21 22 a23 a24; a31 a32 a33 a34; a41 a42 a43 a44];

A=[a11 a12 a13 a14; a21 a22 a23 a24; a31 a32 a33 a34; a41 a42 a43 a44];

B=diag([b11 b22 b33 b44]);

C=diag([c11 c22 c33 c44]);

D=diag([d11 d22 d33 d44]);

%invD=inv(D);

%invA=inv(A);

%B2=[b11 b22 b33; b11 b22 b33; b11 b22 b33];

%inv(B)*B2

E=[A B; C D];

invE=inv(E);