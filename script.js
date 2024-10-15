function isent(form) {
    const gamma = parseFloat(form.g.value); // Get the gamma value (1.4 default)
    const inputType = form.i.value; // Get the input type selected
    const value = parseFloat(form.v.value); // Get the input value

    console.log(gamma);
    console.log(inputType);
    console.log(value);


    if (inputType == "Mach number") {
        let mach = value; 
        isentMach(mach, gamma, form)
    } else if (inputType == "T/T0") {
        let mach = Math.sqrt((2/(gamma - 1))*((1/value) - 1))
        isentMach(mach, gamma, form)
    } else if (inputType == "p/p0") {
        let power = (gamma - 1)/gamma
        let mach = Math.sqrt((2/(gamma - 1)) * ((Math.pow(value, -power)) - 1))
        console.log(mach)
        isentMach(mach, gamma, form)
    } else if (inputType == "rho/rho0") {
        if(value < 0 || value > 1){
            alert("Value should be between 0 and 1")
        }else{
            let power = (gamma - 1)
            let mach = Math.sqrt((2/(gamma - 1)) * ((Math.pow(value, -power)) - 1))
            isentMach(mach, gamma, form)
        }
    } else if (inputType == "A/A* (sub)") {

    }

}

function isentMach(value, gamma, form) {
        let mach, tt0, pp0, rr0, mu, nu, aas;
        mach = value;
        tt0 = 1 / (1 + ((gamma - 1) / 2) * mach * mach);
        pp0 = Math.pow((1 + ((gamma - 1) / 2) * mach * mach), (-(gamma) / (gamma - 1)));
        rr0 = Math.pow((1 + ((gamma - 1) / 2) * mach * mach), (-1 / (gamma - 1)));
        rrs = Math.pow((1 + (gamma - 1) / 2) / (1 + ((gamma - 1) / 2) * mach * mach), 1 / (gamma - 1));
        pps = Math.pow((1 + (gamma - 1) / 2) / (1 + ((gamma - 1) / 2) * mach * mach), gamma / (gamma - 1));
        tts = Math.pow((1 + (gamma - 1) / 2) / (1 + ((gamma - 1) / 2) * mach * mach), 1);
        if (mach > 1) {
            const term1 = Math.sqrt((gamma + 1) / (gamma - 1)) * Math.atan(Math.sqrt((mach * mach - 1) * ((gamma - 1) / (gamma + 1))));
            const term2 = Math.atan(Math.sqrt(mach * mach - 1));
            nu = (term1 - term2) * (180 / Math.PI);
            mu = Math.asin(1 / mach) * (180 / Math.PI);
            const term = (gamma + 1) / 2;
            aas = (1 / mach) * Math.pow(((1 / term) * ((1 + ((gamma - 1) / 2) * mach * mach))), (gamma + 1) / (2 * (gamma - 1))); 
            form.nu.value = nu.toFixed(5);
            form.mu.value = mu.toFixed(5);   
        } else {
            // Set nu and mu back to blank for Mach number <= 1
            nu = "";
            mu = "";
            form.nu.value = nu;
            form.mu.value = mu;
            const term = (gamma + 1) / 2;
            aas = (1 / mach) * Math.pow(((1 / term) * ((1 + ((gamma - 1) / 2) * mach * mach))), (gamma + 1) / (2 * (gamma - 1)));    
        }
        
        form.m.value = mach.toFixed(5);
        form.tt0.value = tt0.toFixed(5);
        form.pp0.value = pp0.toFixed(5);
        form.rr0.value = rr0.toFixed(5);
        form.rrs.value = rrs.toFixed(5);
        form.pps.value = pps.toFixed(5);
        form.tts.value = tts.toFixed(5);
        form.aas.value = aas.toFixed(5);
}


function nsr(form){
    const gamma = parseFloat(form.g.value); // Get the gamma value (1.4 default)
    const inputType = form.i.value; // Get the input type selected
    const value = parseFloat(form.v.value); // Get the input value

    console.log(gamma);
    console.log(inputType);
    console.log(value);

    if(inputType == "M1"){
        let m1 = value;
        nsrMach(m1, form, gamma)
    } else if (inputType == "M2") {
        let m1 = Math.sqrt(((value * value * (gamma - 1)) + 2) / (2 * gamma * value * value - (gamma - 1)));
        nsrMach(m1, form, gamma)
    } else if (inputType == "p2/p1") {
        let m1 = Math.sqrt((value * (gamma + 1) + (gamma - 1)) / (2 * gamma))
        nsrMach(m1, form, gamma)
    } else if (inputType == "T2/T1") {
        let m1 = nsrMachFromTemperature(value, form, gamma)
        nsrMach(m1, form, gamma)
    } else if (inputType == "p02/p01") {
        let m1 = nsrMachFromPressureRatio(value, form, gamma)
        nsrMach(m1, form, gamma)
    }
}

function nsrMach(mach1, form, gamma){
    let m1, m2, p02p01, p02p1, p2p1, r2r1, t2t1;
    m1 = mach1;
    m2 = Math.sqrt((((gamma - 1) * m1 * m1) + 2)/((2 * gamma * m1 * m1) - (gamma - 1)))
    p2p1 = (2 * gamma * m1 * m1 - (gamma - 1))/(gamma + 1)
    t2t1 = (2 * gamma * m1 * m1 - (gamma - 1)) * (2 + (gamma - 1) * m1 * m1) / ((gamma + 1)*(gamma + 1)*m1*m1)
    r2r1 = ((gamma + 1) * m1 * m1) / (2 + (gamma - 1) * m1 * m1)

    let term1 = Math.pow(((gamma + 1)/2) * m1 * m1 / (1 + ((gamma - 1)/2) * m1 *m1), gamma/(gamma - 1))
    let term2 = Math.pow((2 * gamma * m1 * m1 - (gamma - 1))/(gamma + 1), -1/(gamma - 1))
    p02p01 = term1 * term2


    let term3 = Math.pow((gamma + 1) * m1 * m1 / 2, gamma/(gamma - 1))
    let term4 = Math.pow((2 * gamma * m1 * m1 - (gamma - 1)) / (gamma + 1), 1 / (gamma - 1))
    p02p1 = term3 / term4

    form.m1.value = m1.toFixed(5);
    form.m2.value = m2.toFixed(5);
    form.p2p1.value = p2p1.toFixed(5);
    form.t2t1.value = t2t1.toFixed(5);
    form.r2r1.value = r2r1.toFixed(5);
    form.p02p01.value = p02p01.toFixed(5);
    form.p02p1.value = p02p1.toFixed(5);

}

function nsrMachFromTemperature(t2t1, form, gamma) {
    let m1 = 1; // Initial guess
    let epsilon = 0.00001;
    let maxIterations = 1000;
    let iteration = 0;
    let f, f_prime, m1_next;

    function f_m1(m1) {
        return ((2 * gamma * m1 * m1 - (gamma - 1)) * (2 + (gamma - 1) * m1 * m1)) / ((gamma + 1) * (gamma + 1) * m1 * m1) - t2t1;
    }

    // Derivative of the function f(m1) for Newton-Raphson method
    function f_prime_m1(m1) {
        let numerator1 = 4 * gamma * m1 * (2 + (gamma - 1) * m1 * m1);
        let numerator2 = (2 * gamma * m1 * m1 - (gamma - 1)) * (2 * (gamma - 1) * m1);
        let denominator = (gamma + 1) * (gamma + 1) * m1 * m1 * m1;
        return (numerator1 + numerator2) / denominator;
    }

    while (iteration < maxIterations) {
        f = f_m1(m1);
        f_prime = f_prime_m1(m1);
        m1_next = m1 - f / f_prime;

        if (Math.abs(m1_next - m1) < epsilon) {
            break;
        }

        m1 = m1_next;
        iteration++;
    }

    return m1
}


function osr(form) {
    // Input values
    let gamma = parseFloat(form.g.value); // Gamma (specific heat ratio)
    let M1 = parseFloat(form.m.value);   
    let value = parseFloat(form.a.value);
    const inputType = form.i.value;

    if(inputType == "Turn angle (weak shock)"){
        let delta = value
        osrDelta(M1, gamma, delta, form)
    } else if(inputType == "Turn angle (strong shock)") {
        let delta = value
        osrStrongDelta(M1, gamma, delta, form)
    } else if(inputType == "Wave angle"){
        let beta = value
        let betaRad = beta * Math.PI / 180;
    
        let term1 = 2 * (1 / Math.tan(betaRad));
        let term2 = (M1 * M1 * (gamma + Math.cos(2 * betaRad)) + 2); 
        let term3 = (M1 * M1 * Math.sin(betaRad) * Math.sin(betaRad) - 1); 
    
        let tan_delta = term1 * (term3 / term2);
    
        let delta = Math.atan(tan_delta) * 180/ Math.PI;
        console.log(delta)
        osrDelta(M1, gamma, delta, form)
    } else if(inputType == "M1n") {
        let betaRad = Math.asin(value / M1) 

        let term1 = 2 * (1 / Math.tan(betaRad));
        let term2 = (M1 * M1 * (gamma + Math.cos(2 * betaRad)) + 2);
        let term3 = (M1 * M1 * Math.sin(betaRad) * Math.sin(betaRad) - 1);
    
        let tan_delta = term1 * (term3 / term2);
    
        let delta = Math.atan(tan_delta) * 180/ Math.PI;
        console.log(delta)
        osrDelta(M1, gamma, delta, form)
    }
}

function osrDelta(M1, gamma, delta, form){
    delta = delta * Math.PI / 180;     

    function f_beta(beta) {
        return Math.tan(delta) - 2 * (1 / Math.tan(beta)) * ((M1 * M1 * Math.sin(beta) * Math.sin(beta) - 1) / (M1 * M1 * (gamma + Math.cos(2 * beta)) + 2));
    }

    function f_prime_beta(beta) {
        let deltaBeta = 0.0001; 
        return (f_beta(beta + deltaBeta) - f_beta(beta)) / deltaBeta;
    }

    // Newton-Raphson to solve for beta (wave angle)
    let beta = delta; 
    let epsilon = 0.00001;
    let maxIterations = 100;
    let iteration = 0;
    while (iteration < maxIterations) {
        let f = f_beta(beta);
        let f_prime = f_prime_beta(beta);
        let beta_next = beta - f / f_prime;

        if (Math.abs(beta_next - beta) < epsilon) {
            beta = beta_next;
            break;
        }

        beta = beta_next;
        iteration++;
    }


    let beta_deg = beta * 180 / Math.PI;
    
    let M1n = M1 * Math.sin(beta);

    let M2n = Math.sqrt(((gamma - 1) * M1n * M1n + 2) / (2 * gamma * M1n * M1n - (gamma - 1)));
    let M2 = M2n / Math.sin(beta - delta); // Oblique shock relation for M2

    let p2p1 = 1 + (2 * gamma / (gamma + 1)) * (M1n * M1n - 1);
    let r2r1 = ((gamma + 1) * M1n * M1n) / ((gamma - 1) * M1n * M1n + 2);
    let t2t1 = p2p1 / r2r1;

    let term1 = Math.pow(((gamma + 1) / 2) * M1n * M1n / (1 + ((gamma - 1) / 2) * M1n * M1n), gamma / (gamma - 1));
    let term2 = Math.pow((2 * gamma * M1n * M1n - (gamma - 1)) / (gamma + 1), -1 / (gamma - 1));
    let p02p01 = term1 * term2;

    form.m2.value = M2.toFixed(5);
    form.delta.value = (delta * 180 / Math.PI).toFixed(5); 
    form.beta.value = beta_deg.toFixed(5);                 
    form.p2p1.value = p2p1.toFixed(5);                     
    form.r2r1.value = r2r1.toFixed(5);                     
    form.t2t1.value = t2t1.toFixed(5);                     
    form.p02p01.value = p02p01.toFixed(5);                 
    form.m1n.value = M1n.toFixed(5);                       
    form.m2n.value = M2n.toFixed(5);                       
}

function osrStrongDelta(M1, gamma, delta, form) {
    delta = delta * Math.PI / 180;

    // Find the strong shock wave angle beta using Newton-Raphson method
    function f_beta(beta) {
        return Math.tan(delta) - 2 * (1 / Math.tan(beta)) * ((M1 * M1 * Math.sin(beta) * Math.sin(beta) - 1) / (M1 * M1 * (gamma + Math.cos(2 * beta)) + 2));
    }

    function f_prime_beta(beta) {
        let delta_beta = 0.0001;
        return (f_beta(beta + delta_beta) - f_beta(beta)) / delta_beta;
    }

    // Newton-Raphson to solve for beta (wave angle)
    let beta = Math.PI / 2;
    let epsilon = 0.00001;
    let maxIterations = 100;
    let iteration = 0;
    while (iteration < maxIterations) {
        let f = f_beta(beta);
        let f_prime = f_prime_beta(beta);
        let beta_next = beta - f / f_prime;

        if (Math.abs(beta_next - beta) < epsilon) {
            beta = beta_next;
            break;
        }

        beta = beta_next;
        iteration++;
    }

    let beta_deg = beta * 180 / Math.PI;

    let M1n = M1 * Math.sin(beta);

    let M2n = Math.sqrt(((gamma - 1) * M1n * M1n + 2) / (2 * gamma * M1n * M1n - (gamma - 1)));
    let M2 = M2n / Math.sin(beta - delta); 

    let p2p1 = 1 + (2 * gamma / (gamma + 1)) * (M1n * M1n - 1);

    let r2r1 = ((gamma + 1) * M1n * M1n) / ((gamma - 1) * M1n * M1n + 2);

    let t2t1 = p2p1 / r2r1;

    let term1 = Math.pow(((gamma + 1) / 2) * M1n * M1n / (1 + ((gamma - 1) / 2) * M1n * M1n), gamma / (gamma - 1));
    let term2 = Math.pow((2 * gamma * M1n * M1n - (gamma - 1)) / (gamma + 1), -1 / (gamma - 1));
    let p02p01 = term1 * term2;

    form.m2.value = M2.toFixed(5);
    form.delta.value = (delta * 180 / Math.PI).toFixed(5);
    form.beta.value = beta_deg.toFixed(5);                
    form.p2p1.value = p2p1.toFixed(5);                    
    form.r2r1.value = r2r1.toFixed(5);                    
    form.t2t1.value = t2t1.toFixed(5);                    
    form.p02p01.value = p02p01.toFixed(5);               
    form.m1n.value = M1n.toFixed(5);                      
    form.m2n.value = M2n.toFixed(5);                      
}
