% Electrostatics
onesix = 0.16666666667;
twothirds = 0.6666666667;
AKPERP = 2.405;
K8 = 8192;
NGMAX = K8;
NXV = K8;
NVY = NXV / 2;
NSPM = 8;
HISTMAX = 512;
NVBINMAX = 1024;
MAXPARTICLES = 920001;
m_ng;
m_dx;
m_x;
m_l, m_dt, m_la, m_ael, m_epsi, m_rho0, m_a1, m_a2, m_e0, m_w0, m_ese_hist, m_p, m_ecconst;
m_kes_hist, m_pxs_hist, m_esem_hist, m_ms, m_qs, m_ts, m_nms;
m_ese, m_kes, m_pxs, m_esem, m_ke, m_te, m_vx, m_vy;
m_x_array, m_t_array, m_k_array, m_rho, m_phi, m_phik, m_e;
m_acc, m_dvbin, m_v_array, m_vbint, m_vbin, m_vbinstart, m_vbin_inst;
m_t;
m_it, m_nt, m_ith, m_ithl, m_hist_hi, m_interval, m_nsp, m_localnsp, m_accum, m_mmax, m_ng1, m_iw, m_ec, m_k_hi;
m_ins = [1:NSPM + 2];
m_np = [1:NSPM + 2];
m_vbins = [1:NSPM + 2];
m_nvbin = [1:NSPM];

function retval = init(il1, il2, m, q, t, nm, vbin1, vbin2, dvb, vstart, nvbin)
    retval = true;
    wp = 1.0, wc = 0.0, qm = -1.0, vt1 = 0.0, vt2 = 0.0, v0 =0.0, x1 = 0.0, v1 = 0.0, thetax = 0.0, thetav = 0.0;
    ddx, x0, vmax, dv, vvnv2, vv, fv, df, xs, xsi, theta, lg, vupper, vlower;
    ngr, j, i1, i2;
    n = 128, nlg = 1, nv2 = 0, mode = 1, nbins = 100;
    a_char = [1:80];
    
    dvb = 1.0;

    %{
    while (fscanf(InputDeck, "%d %d %d %d ", &n, &nv2, &nlg, &mode) < 4)  // added nbins
    fscanf(InputDeck,"%s", a_char);
    while (fscanf(InputDeck, "%g %g %g %g %g %g", &wp, &wc, &qm, &vt1, &vt2, &v0) < 6)
    fscanf(InputDeck,"%s", a_char);
    while (fscanf(InputDeck,"%g %g %g %g", &x1, &v1, &thetax, &thetav) < 4)
    fscanf(InputDeck,"%s", a_char);
    while(fscanf(InputDeck," %d %g %g ",&nbins,&vlower,&vupper) < 3)
    fscanf(InputDeck,"%s", a_char);
    %}
    fprintf(' n = %4d  nv2 = %4d  nlg = %4d  mode = %4d \n', n, nv2 ,nlg, mode);
    fprintf(' wp = %6.3f   wc = %6.3f   qm = %6.3f \n', wp, wc, qm);
    fprintf(' vt1 = %6.3f  vt2 = %6.3f  v0 = %6.3f \n', vt1, vt2, v0);
    fprintf(' x1 = %6.3f   v1 = %6.3f   thetax = %6.3f   thetav = %6.3f \n', x1, v1, thetax, thetav);
    fprintf(' nbins = %4d   vlower = %6.3f  vupper = %6.3f ', nbins, vlower, vupper);

    t = tan(-0.5 * wc * m_dt);
    il2 = il1 + n;

    q = m_l * wp * wp / (m_epsi * n * qm);
    m = q / qm;
    nm = n * m;
    ngr = n / nlg;
    lg = m_l / nlg;
    ddx = m_l / n;


    if (nbins > NVBINMAX)
        nbins = NVBINMAX;
    end
    if (nbins < 2)
        nbins = 2;
    end
    nvbin = nbins;
    vbin2 = vbin1 + nbins;
    if (vupper - vlower < 0.0)
        fprintf('\nError in INIT: vupper must be > vlower!');
        retval = false;
        return;
    end
    if (vt1 < 0 || vt2 < 0)
        fprintf('\nError in INIT: can't have negative thermal voltages!");
        retval = false;
        return;
    end
    if (vupper-vlower > 0.0)
        vstart = vlower;
        dvb = (vupper - vlower) / nbins;
    else
        if (vt1 + vt2 > 0.0)
            vstart = v0 - 5.0 * (vt1 + vt2);
            dvb = 10.0 * (vt1 + vt2) / nbins;  % so that the distribution goes from v0-5*vt to v0+5*vt
        elseif (fabs(v0) > 0)
            if (v0 > 0)
                vstart = 0.0;
                dvb = 2 * v0 / nbins;
            else
                vstart = 2.0 * v0;
                dvb = -2 * v0 / nbins;
            end
        else
            vstart = 0.0;
            dvb = 1 / nbins;
        end
    end

    % setup v_array for this species
    for i = vbin1:vbin2
        m_v_array[i] = (vstart + (i - vbin1 + 0.5) * dvb);
    end

    for i = 1:ngr
        i1 = i - 1 + il1;
        x0 = (i - 0.5) * ddx;
        m_x[i1] = x0;
        m_vx[i1] = v0;
    end
    if (vt2 ~= 0.0)
        vmax = 5.0 * vt2;
        dv = 2.0 * vmax / (n - 1);
        vvnv2 = 1.0;
        m_x[il1] = 0.0;
        for i = 2:n
            vv = ((i - 1.5) * dv - vmax) / vt2;
            if (nv2 ~= 0)
                vvnv2 = vv ^ nv2;
            end
            fv = vvnv2 * exp(-0.5 * vv * vv);
            i1 = i - 1 + il1;
            m_x[i1] = m_x[i1 - 1] + ((fv >= 0.0) ? fv : 0.0);
            end
        df = m_x[i1]/ngr;
        i1 = il1;
        j = il1;
        for i = 1:ngr
            fv = (i - 0.5) * df;
            while (fv >= m_x[j + 1])
                j++;
                if (j > (il2 - 2))
                    fprintf('distribution function error');
                    retval = false;
                    return
                end
            end
            vv = dv * (j - il1 + (fv - m_x[j]) / (m_x[j + 1] - m_x[j])) - vmax;
            m_vx[i1] += vv;
            i1++;
        end
        xs = 0.0;
        for i = 1:ngr
            i1 = i - 1 + il1;
            m_x[i1] = xs * lg + 0.5 * ddx;
            xsi = 1.0;
            xsi *= 0.5;
            xs -= xsi;
            while (xs >= 0.0)
                xsi *= 0.5;
                xs -= xsi;
            end
            xs += 2.0 * xsi;
        end
        i1 = ngr + il1 - 1;
    end
    if (wc ~= 0.0)
        for i = 1:ngr
            i1 = i - 1 + il1;
            vv = m_vx[i1];
            theta = 2 * PhysConsts::PI * frand();
            m_vx[i1] = vv * cos(theta);
            m_vy[i1] = vv * sin(theta);
        end
    end
    if (nlg ~= 1)
        j = ngr + 1;
        xs = 0.0;
        for i = j:ngr:n
            xs += lg;
            for j = 1:ngr
                i1 = j - 1 + il1;
                i2 = i1 + i - 1;
                m_x[i2] = m_x[i1] + xs;
                m_vx[i2] = m_vx[i1];
                if (wc ~= 0.0)
                    m_vy[i2] = m_vy[i1];
                end
            end
        end
    end
    if (vt1 ~= 0.0)
        for i = 1:1:n
            i1 = i - 1 + il1;
            for j = 0:12
                if (wc ~= 0.0)
                    m_vy[i1] += vt1 * (frand() - 0.5);
                end
                m_vx[i1] += vt1 * (frand() - 0.5);
            end
        end
    end
    for i = 1;1:n
        i1 = i - 1 + il1;
        theta = 2 * pi * mode * m_x[i1] / m_l;
        m_x[i1] += x1 * cos(theta + thetax);
        m_vx[i1] += v1 * sin(theta + thetav);
    end
    setrho(il1, il2 - 1, q, q * n / m_l);
end


function main(argc, argv)
    char a_char[80];
    bool WasInputFileGiven = false;
    QString filename;

    % InputDeck = (!WasInputFileGiven) ? fopen("es1data","r") :  fopen(theInputFile,"r");
    % if (!InputDeck) {
    if (0)
        fprintf('Can''t find input file %s', argv[1]);
        fprintf('Correct syntax is: ES1 -i file.inp');
    else
%{
        read lines until we get to numbers
        while (fscanf(InputDeck,"%d %g %g %d %d %g %d", &nsp, &l, &dt, &nt, &mmax, &la, &accum) < 7)
            fscanf(InputDeck, "%s", a_char);
        // note: la is l/a

        while (fscanf(InputDeck," %d %d %d %g %g %g %g %g", &ng, &iw, &ec, &epsi, &a1, &a2, &e0, &w0) < 8)
            fscanf(InputDeck, "%s", a_char);
        %}
        if (m_nsp > NSPM)
            fprintf('Number of species nsp cannot exceed NSPM');
            return;
        end

        if (m_accum < 0)
            fprintf('Error:  accum can't be negative!');
            return;
        end

        if (m_ec ~= 1 || m_ec ~= 0)
            fprintf('Error:  What are you thinking?  There are only two possible values of ec');
            fprintf('0 and 1.  %d is not 0 or 1.', m_ec);
            return;
        end
        if (not m_iw && m_ec)
            fprintf('Error:  There IS no energy-conserving algorithm for NGP');
            return;
        end

        m_ecconst = (m_ec) ? 0.5 : 0.0;
        if (m_iw > 3 || m_iw < 0)
            fprintf('Error:  bad iw flag!  Please check your input deck!');
            return;
        end

        if (m_ng > NGMAX) {
            fprintf('Number of grids ng cannot exceed NGMAX');
            return;
        end
        m_dx = m_l / m_ng;
        m_ng1 = m_ng + 1;
        m_k_hi = m_ng / 2;

        % Allocating space for arrays
        m_nms = zeros(m_nsp + 1);
        m_ms = zeros(m_nsp + 1];
        m_qs = zeros(m_nsp + 1];
        m_ts = zeros(m_nsp + 1];

        m_x_array = zeros(m_ng + 1];
        for i = 0:m_ng
            m_x_array[i] = i * m_dx;
        end
        m_rho = zeros(m_ng + 3];
        m_phi = zeros(m_ng + 2];
        m_phik = zeros(m_ng + 2];
        m_k_array = zeros(m_ng];
        for i = 0:m_k_hi
            m_k_array[i] = i * 2 * pi / m_l;
        end
        m_e = zeros(m_ng + 2];
        m_acc = zeros(m_ng + 3];

        m_t_array = zeros(HISTMAX];
        m_ese = zeros(HISTMAX];
        m_ke = zeros(HISTMAX];
        m_te = zeros(HISTMAX];

        m_kes_hist = zeros(m_nsp + 1];
        m_pxs_hist = zeros(m_nsp + 1];
        m_esem_hist = zeros(m_mmax + 1];

        m_kes = zeros([m_nsp];
        for i = 0:m_nsp
            m_kes[i] = zeros(HISTMAX];
        end
        m_pxs = zeros(m_nsp];
        for i = 0:m_nsp
            m_pxs[i] = zeros(HISTMAX];
        end
        m_esem = zeros(m_mmax];
        for i = 0:m_mmax
            m_esem[i] = zeros(HISTMAX];
        end
        
        m_x = zeros(MAXPARTICLES];
        m_vx = zeros(MAXPARTICLES];
        m_vy = zeros(MAXPARTICLES];

        m_vbint = zeros(NVBINMAX];
        m_vbin = zeros(m_nsp * NVBINMAX];
        m_vbin_inst = zeros(m_nsp * NVBINMAX];
        for i = 0:m_nsp * NVBINMAX
            m_vbin_inst[i]= 0.0;
        end
        m_dvbin = zeros(m_nsp + 1];
        m_vbinstart = zeros(m_nsp + 1];
        m_v_array = zeros(NVBINMAX * m_nsp];

        if (!m_x || !m_vx || !m_vy || !m_vbint || !m_vbin_inst || !m_dvbin || !m_v_array )
            fprintf('START: Could not get enough memory for x or v''s.');
            return;
        end

        fprintf('nsp = %2d     l = %8.5f', m_nsp, m_l);
        fprintf('dt = %4.2f    nt = %4d', m_dt, m_nt);
        fprintf('ng = %5d   iw = %2d   ec = %2d  accum = %4d', m_ng, m_iw, m_ec, m_accum);
        fprintf('epsi = %4.2f  a1 = %4.2f  a2 = %4.2f', m_epsi, m_a1, m_a2);

        for i = 1:m_nsp
            if (!init(&m_ins[i], &m_ins[i + 1], &m_ms[i], &m_qs[i], &m_ts[i], &m_nms[i], &m_vbins[i], &m_vbins[i + 1], &m_dvbin[i], &m_vbinstart[i], &m_nvbin[i]))
                return;
            end
        end

        % added vbins to param list
        % fclose(InputDeck);

        for i = 1:m_nsp
            m_np[i] = m_ins[i + 1] - m_ins[i];
        end
        for i = 0:m_nsp
            for j = 0:HISTMAX
                m_kes[i][j] = 0.0;
                m_pxs[i][j] = 0.0;
            end
        end
        m_rho[1] += m_rho[m_ng1];   %These resolve the periodic boundary conditions.
        m_rho[2] += m_rho[m_ng + 2];
        m_rho[m_ng] += m_rho[0];
        m_rho[0] = m_rho[m_ng];
        m_rho[m_ng + 2] = m_rho[2];
        m_rho[m_ng1] = m_rho[1];
        fields(0);

        for i = 1:m_nsp
            setv(m_ins[i], m_ins[i + 1] - 1, m_qs[i], m_ms[i], m_ts[i], m_pxs_hist[i], m_kes_hist[i]);
        end
        % scale all the velocities  properly
        m_dvbin[i] *= m_dt / m_dx;
        m_vbinstart[i] *= m_dt / m_dx;
        startvel();
    end
end
 

function [retval] = sign(a, b)
    if (b >= 0.0)
        if (a >= 0.0)
            retval = a;
        else
            retval = -a;
        end
    else
        if (a < 0.0)
            retval = a;
        else
            retval = -a;
        end        
    end
end
%{
void startvel();
void accel(int, int, double, double, double, double *, double *);
void cpft(double [], double [], int, int, double);
void rpft2(double [], double [] ,const int, const int);
void rpfti2(double [], double [],const int, const int);
void fields(const int);
void move(const int ilp, const int, const double);
void setv(const int, const int, const double, const double, const double, double *, double *);
void setrho(const int, const int, const double, const double);
double frand();
void velocity();
void history();

double sqr(const double val) { return val * val; }
double cube(const double val) { return val * val * val; }

// member properties
int m_ng;
double m_dx;
double *m_x;
double m_l, m_dt, m_la, m_ael, m_epsi, m_rho0, m_a1, m_a2, m_e0, m_w0, m_ese_hist, m_p, m_ecconst;
double *m_kes_hist, *m_pxs_hist, *m_esem_hist, *m_ms, *m_qs, *m_ts, *m_nms;
double *m_ese, **m_kes, **m_pxs, **m_esem, *m_ke, *m_te, *m_vx, *m_vy;
double *m_x_array, *m_t_array, *m_k_array, *m_rho, *m_phi, *m_phik, *m_e;
double *m_acc, *m_dvbin, *m_v_array, *m_vbint, *m_vbin, *m_vbinstart, *m_vbin_inst;
double m_t;
int m_it, m_nt, m_ith, m_ithl, m_hist_hi, m_interval, m_nsp, m_localnsp, m_accum, m_mmax, m_ng1, m_iw, m_ec, m_k_hi;
int m_ins[NSPM + 2], m_np[NSPM + 2], m_vbins[NSPM + 2], m_nvbin[NSPM];
}%