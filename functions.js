
function factorial(n) {
	++n;
	var r = 1;
  while(n-->1)
  	r*=n;
   return r;
};

function PI(n) {
	n = n || 30;
    var k = Math.sqrt(2)*2/9801;
    var sum = 0;
    for(var i=0; i<n; i++) {
          sum += (factorial(4*i)*(1103+26390*i))/
      (Math.pow(factorial(i), 4)*Math.pow(396, 4*i));
    }
    return 1/(sum*k);
};

function qlog(x, e=1e-16) {
    var a, b, c, d, de, ans, prev, k;
    a = x-1;
    b = x+1;
    c = (a*a)/(b*b);
    d = 2*a/b;
    ans = 0;
    k = 0;
    prev = 0;
    do {
        ans += (1/(2*k+1))*Math.pow(c, k);
        de = ans - prev;
        prev = ans;
        k++;
    }
    while(de >= e)
    return ans*d;
};

//arithmatic mean
function M(x, y) {
    var a, b, d, e=1e-20;
    do {
        a = (x+y)/2;
        b = Math.sqrt(x*y);
        d = a-b;
        x = a;
        y = b;
    }
    while(d >= e)
    return x;
};

function calc_m(p, x) {
    var m = 1;
    var p2 = p/2;
    while(x*Math.pow(2, m) < Math.pow(2, p2)) 
        m++;

    return m;
};

function log(x, p=25) {
    var m, LN2, pi;
    LN2 = qlog(2);
    m = calc_m(p, x);
    pi = PI();
    return pi/(2*M(1, Math.pow(2, 2-m)/x))-m*LN2;
};

function nthroot(x, n) {
    var xk, xk1, d, e;
    e = 1e-16;
    xk = 3*x/(2*n);
    do {
        xk1 = 1/n*((n-1)*xk+x/pow(xk, n-1));
        d = xk-xk1;
        xk = xk1;
    }
    while(d > e)
    return xk;
};

//https://en.wikipedia.org/wiki/Exponentiation_by_squaring
function pow(x, n) {
    if(n > 0 && n < 1) {
        //convert to a quick fraction
        var f = fraction(n);
        return pow(nthroot(x, f[1]), f[0]);
    }
    if(n < 0) {
        x = 1/x;
        n = -n;
    }
    if(n === 0)
        return 1;
    var y = 1;
    while(n > 1) {
        if(n%2 === 0) {
            x = x*x;
            n = n/2;
        }
        else {
            y = x*y;
            x = x*x;
            n = (n-1)/2;
        }
    }
    return x*y;
};

function abs(x) {
    if(x >= 0)
        return x;
    return -x;
}

function gcd(a,b) {
    a = abs(a);
    b = abs(b);
    if (b > a) 
        return gcd(b, a);
    while (true) {
        if (b === 0) 
            return a;
        a %= b;
        if (a === 0) 
            return b;
        b %= a;
    }
}

function qfrac(n) {
    var a, q;
    a = 1;
    while(n%1 !== 0) {
        n *= 10;
        a *= 10;
    }
    q = gcd(a, n);
    return [n/q, a/q];
}

function sign(x) {
    return x/abs(x);
}

function fraction(x, e=0.0019){
    var w, f, n, s, stack, dec, num;
    s = sign(x);
    stack = 0;
    w = Math.floor(x);
    n = x-w;
    /*recursive function that transforms the fraction*/
    function convert(x){
        stack++;
        var intgr = Math.floor(x); //get the integer part of the number
        dec = (x - intgr); //get the decimal part of the number
        if(dec < e || stack > 20) 
            return [intgr,1]; //return the last integer you divided by
        num = convert(1/dec); //call the function again with the inverted decimal part
        return[intgr*num[0]+num[1],num[0]];
    }
    f = convert(n); 
    f[0] += w*f[1];
    f[0] *= s;
    return f;
}
