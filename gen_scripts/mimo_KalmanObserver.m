clc;
close all;
clear variables;
format longEng

%%

A = [ ...
   0,  1,  2; ...
  -2, -3,  0; ...
  -2,  0, -1; ...
  ];
B = [ ...
  2, 1;...
  0, 2;...
  1, 3;...
  ];
C = [ ...
  1,  2, -1; ...
  2, -1,  3; ...
  ];


sys_ref = c2d(ss(A,B,C,0),0.2);

%% Kalman

sys = sys_ref;

A_kal = [ ...
      sys.A,   sys.B; ...
 zeros(2,3),  eye(2); ...
  ];

B_kal = [ ...
       sys.B; ...
  zeros(2,2); ...
  ];

C_kal = [ sys.C, zeros(2,2) ];

nn_kal = size(A_kal,2);
ni_kal = size(B_kal,2);
no_kal = size(C_kal,1);

C_int = C_kal(:,:);


%%
R_kal = blkdiag(1E-2, 1E-2);
Q_kal = blkdiag(B_kal(1:3,:)*B_kal(1:3,:)'* 1E-3,  1E-16,  1E-16);

P0    = eye(nn_kal);
X0    = zeros(nn_kal,1);

%%  Regler

  A_reg = [ ...
   sys_ref.A, zeros(3,2); ...
  -sys_ref.C,     eye(2); ...
  ];

B_reg = [ ...
  sys_ref.B; ...
 zeros(2,2); ...
  ];

nn = 3;
ni = 2;
no = 2;

%%
Q_reg = blkdiag(1,1,1,1,1)*1E0;
R_reg = blkdiag(1,1)*1E3;
[hT_star, ~, e] = dlqr(A_reg, B_reg, Q_reg, R_reg, zeros(nn+ni,no));

e

%%

hT_kal = [hT_star(:, 1:nn), zeros(2,2)]
hT_int = hT_star(:, nn+1:end)

A_int = eye(no,no);
B_int = [ ...
  -hT_int^-1, hT_int^-1, eye(ni); ...
  ];


%%
Ts = 1/ 1E4;
N = 2E3;

t = (0:N-1)' * Ts;
r = zeros(2, numel(t));

r(1,0.05<t) = 1;
r(2,0.05<t) = 1;

r(1,0.10<t) = 1;
r(2,0.10<t) = 0;

r(1,0.15<t) = 0;
r(2,0.15<t) = 1;

x      = nan(nn_kal, N);
x_int  = nan(     2, N);
y      = nan(no_kal, N);

u_star = nan(2, N);
u_dash = nan(2, N);
u      = nan(2, N);
e      = nan(2, N);

y_ref  = nan(no_kal, N);
y_meas = nan(no_kal, N);
x_sim  = nan(size(sys_ref.A,1), N);
y_sim  = nan(size(sys_ref.C,1), N);

y_noise = randn(size(sys_ref.C,1), N) * 1E-3;


A_sim = sys_ref.A;
B_sim = sys_ref.B;
C_sim = sys_ref.C;

PP     = P0;
x_kal_ = X0;
x_int_ = zeros(2,1);
x_sim_ = zeros(size(A_sim,1), 1);

for idx = 1:numel(t)
    K = PP * C_kal' * (C_kal * PP * C_kal' + R_kal)^-1;
    P_ = (eye(nn_kal) - K*C_kal) * PP;
    PP = A_kal * P_ * A_kal' + Q_kal;
    
    y_meas_ = C_sim * x_sim_;
    y_ref(:, idx) = y_meas_;
    y_meas_ = y_meas_ + y_noise(:,idx);
    y_meas(:, idx) = y_meas_;
    
    y_ = C_kal*x_kal_;
    x_kal_ = x_kal_ + K * (y_meas_ - y_);
    
    y(:, idx)    = C_kal*x_kal_;
    x(:, idx)    = x_kal_;
    x_int(:,idx) = x_int_;

    u_star_ = -hT_kal*x_kal_ - hT_int*x_int_;
    u_star(:, idx) = u_star_;
    
    u_dash_ = clip(u_star_, [1;1]*-2.0, [1;1]*2.0);
    u_dash(:, idx) = u_dash_;
    
    u_ = u_dash_; round(u_dash_, 3);

    u(:, idx) = u_;
    
    r_ = r(:,idx);
    e_ = r_ - C_int * x_kal_;
    e(:,idx) = e_;
    
    x_kal_ = A_kal * x_kal_ + B_kal*u_;
    x_int_ = A_int * x_int_ + B_int*[u_dash_; u_star_; e_];
    
    x_sim_ = A_sim * x_sim_ + B_sim * u_;
end

%%
figure(); hold on; grid on;
    pl_y = subplot(311); hold on; grid on;
    pl_u = subplot(312); hold on; grid on;
    pl_n = subplot(313); hold on; grid on;

  
    plot(pl_y, t, y_meas(1,:), 'r.', LineWidth=1);
    plot(pl_y, t, y(1,:),      'r-', LineWidth=1);

    plot(pl_y, t, y_meas(2,:), 'b.', LineWidth=1);
    plot(pl_y, t, y(2,:),      'b-', LineWidth=1);
    
    plot(pl_y, t, r, 'k', LineWidth=1);
    
    plot(pl_u, t, u_star(:,:), 'b.', LineWidth=1);
    plot(pl_u, t, u_dash(:,:), 'r-', LineWidth=1);
    plot(pl_u, t, u(:,:),      'g-', LineWidth=1);
    
    plot(pl_n, t, x(end-1:end,:), 'm-', LineWidth=1);

%%

  sel = t < 0.5;
  test_data.y_meas = y_meas(:,sel);
  test_data.u_star = u_star(:,sel);
  test_data.u_dash = u_dash(:,sel);
  test_data.u      = u(:,sel);
  test_data.y      = y(:,sel);
  test_data.x      = x(:,sel);
  test_data.r      = r(:,sel);

%%
  export_file_name   = "tst_data_calc_MIMO_KalmanObserver";
  namespace_location = "test_data::parameter";

  str_hpp = "";
  str_hpp = str_hpp + sprintf("#pragma once\n");
  str_hpp = str_hpp + sprintf("#ifndef %s_HPP_INCLUDED\n", strrep(upper(export_file_name), ".", "_"));
  str_hpp = str_hpp + sprintf("#define %s_HPP_INCLUDED\n\n", strrep(upper(export_file_name), ".", "_"));
  
  str_hpp = str_hpp + sprintf("#include <exmath.hpp>\n\n");
  str_hpp = str_hpp + sprintf("#include <controller.hpp>\n\n");
  
  str_hpp = str_hpp + sprintf("namespace %s\n{\n", namespace_location);
  str_hpp = str_hpp + sprintf("  using parameter_t = controller::calculator::MIMO_KalmanObserver<float, 5, 2, 2>::parameter_t;\n");
  str_hpp = str_hpp + sprintf("  extern parameter_t const param;\n");
  if exist("test_data", "var")
    str_hpp = str_hpp + sprintf("  extern exmath::matrix_t<float, %d, %d> const y_meas;\n", size(test_data.y_meas));
    str_hpp = str_hpp + sprintf("  extern exmath::matrix_t<float, %d, %d> const y;\n", size(test_data.y));
    str_hpp = str_hpp + sprintf("  extern exmath::matrix_t<float, %d, %d> const u;\n", size(test_data.u));
    str_hpp = str_hpp + sprintf("  extern exmath::matrix_t<float, %d, %d> const x;\n", size(test_data.x));
  end
  str_hpp = str_hpp + sprintf("} // namespace %s\n\n", namespace_location);
  str_hpp = str_hpp + sprintf("#endif\n\n");

  fid = fopen(sprintf("out/%s.hpp", export_file_name), "w");
    assert(fid);
    fprintf(fid,"%s", str_hpp);
  fclose(fid);

  str = "";
  str = str + sprintf('#include "%s.hpp\"\n\n', export_file_name);
  str = str + sprintf("namespace %s\n{\n", namespace_location);
  str = str + sprintf("  parameter_t const param{\n");
  str = str + sprintf("    .sys{\n");
  str = str + sprintf("      .A%s,\n",  print_matrix(A_kal, "float"));
  str = str + sprintf("      .B%s,\n",  print_matrix(B_kal, "float"));
  str = str + sprintf("      .C%s,\n",  print_matrix(C_kal, "float"));
  str = str + sprintf("      .X0%s\n", print_matrix(X0,    "float"));
  str = str + sprintf("    },\n");
  str = str + sprintf("    .Q%s,\n",    print_matrix(Q_kal, "float"));
  str = str + sprintf("    .R%s,\n",    print_matrix(R_kal, "float"));
  str = str + sprintf("    .P0%s\n",   print_matrix(P0,    "float"));
  str = str + sprintf("  };\n");
  if exist("test_data", "var")
    str = str + sprintf("  exmath::matrix_t<float, %d, %d> const y_meas%s;\n", size(test_data.y_meas), print_matrix(test_data.y_meas, "float"));
    str = str + sprintf("  exmath::matrix_t<float, %d, %d> const y%s;\n", size(test_data.y), print_matrix(test_data.y, "float"));
    str = str + sprintf("  exmath::matrix_t<float, %d, %d> const u%s;\n", size(test_data.u), print_matrix(test_data.u, "float"));
    str = str + sprintf("  exmath::matrix_t<float, %d, %d> const x%s;\n", size(test_data.x), print_matrix(test_data.x, "float"));
  end
  str = str + sprintf("} // namespace %s\n\n", namespace_location);

  fid = fopen(sprintf("out/%s.cpp", export_file_name), "w");
    assert(fid);
    fprintf(fid,"%s", str);
  fclose(fid);

%%
  export_file_name   = "tst_data_calc_MIMO_KalmanIntegralController";
  namespace_location = "test_data::parameter";

  str_hpp = "";
  str_hpp = str_hpp + sprintf("#pragma once\n");
  str_hpp = str_hpp + sprintf("#ifndef %s_HPP_INCLUDED\n", strrep(upper(export_file_name), ".", "_"));
  str_hpp = str_hpp + sprintf("#define %s_HPP_INCLUDED\n\n", strrep(upper(export_file_name), ".", "_"));
  
  str_hpp = str_hpp + sprintf("#include <exmath.hpp>\n\n");
  str_hpp = str_hpp + sprintf("#include <controller.hpp>\n\n");
  
  str_hpp = str_hpp + sprintf("namespace %s\n{\n", namespace_location);
  str_hpp = str_hpp + sprintf("  using parameter_t = controller::calculator::MIMO_KalmanIntegralController<float, 5, 2, 2, 2>::parameter_t;\n");
  str_hpp = str_hpp + sprintf("  extern parameter_t const param;\n");
  if exist("test_data", "var")
    str_hpp = str_hpp + sprintf("  extern exmath::matrix_t<float, %d, %d> const y_meas;\n", size(test_data.y_meas));
    str_hpp = str_hpp + sprintf("  extern exmath::matrix_t<float, %d, %d> const u_star;\n", size(test_data.u_star));
    str_hpp = str_hpp + sprintf("  extern exmath::matrix_t<float, %d, %d> const u_dash;\n", size(test_data.u_dash));
    str_hpp = str_hpp + sprintf("  extern exmath::matrix_t<float, %d, %d> const u;\n", size(test_data.u));
    str_hpp = str_hpp + sprintf("  extern exmath::matrix_t<float, %d, %d> const x;\n", size(test_data.x));
    str_hpp = str_hpp + sprintf("  extern exmath::matrix_t<float, %d, %d> const r;\n", size(test_data.r));
  end
  str_hpp = str_hpp + sprintf("} // namespace %s\n\n", namespace_location);
  str_hpp = str_hpp + sprintf("#endif\n\n");

  fid = fopen(sprintf("out/%s.hpp", export_file_name), "w");
    assert(fid);
    fprintf(fid,"%s", str_hpp);
  fclose(fid);

  str = "";
  str = str + sprintf('#include "%s.hpp\"\n\n', export_file_name);
  str = str + sprintf("namespace %s\n{\n", namespace_location);
  str = str + sprintf("  parameter_t const param{ // ... \n");
  str = str + sprintf("    .observer{\n");
  str = str + sprintf("      .sys{\n");
  str = str + sprintf("        .A%s,\n",print_matrix(A_kal, "float"));
  str = str + sprintf("        .B%s,\n",print_matrix(B_kal, "float"));
  str = str + sprintf("        .C%s,\n",print_matrix(C_kal, "float"));
  str = str + sprintf("        .X0%s\n", print_matrix(X0,    "float"));
  str = str + sprintf("      },\n");
  str = str + sprintf("      .Q%s,\n",print_matrix(Q_kal, "float"));
  str = str + sprintf("      .R%s,\n",print_matrix(R_kal, "float"));
  str = str + sprintf("      .P0%s\n",   print_matrix(P0,    "float"));
  str = str + sprintf("    },\n");
  str = str + sprintf("    .C_int%s,\n",print_matrix(C_int, "float"));
  str = str + sprintf("    .hT_kal%s,\n",print_matrix(hT_kal, "float"));
  str = str + sprintf("    .A_int%s,\n",print_matrix(A_int, "float"));
  str = str + sprintf("    .B_int%s,\n",print_matrix(B_int, "float"));
  str = str + sprintf("    .hT_int%s,\n",print_matrix(hT_int, "float"));
  str = str + sprintf("  };\n");
  if exist("test_data", "var")
    str = str + sprintf("  exmath::matrix_t<float, %d, %d> const y_meas = %s;\n", size(test_data.y_meas), print_matrix(test_data.y_meas, "float"));
    str = str + sprintf("  exmath::matrix_t<float, %d, %d> const u_star = %s;\n", size(test_data.u_star), print_matrix(test_data.u_star, "float"));
    str = str + sprintf("  exmath::matrix_t<float, %d, %d> const u_dash = %s;\n", size(test_data.u_dash), print_matrix(test_data.u_dash, "float"));
    str = str + sprintf("  exmath::matrix_t<float, %d, %d> const u      = %s;\n", size(test_data.u), print_matrix(test_data.u, "float"));
    str = str + sprintf("  exmath::matrix_t<float, %d, %d> const x      = %s;\n", size(test_data.x), print_matrix(test_data.x, "float"));
    str = str + sprintf("  exmath::matrix_t<float, %d, %d> const r      = %s;\n", size(test_data.r), print_matrix(test_data.r, "float"));
  end
  str = str + sprintf("} // namespace %s\n\n", namespace_location);

  fid = fopen(sprintf("out/%s.cpp", export_file_name), "w");
    assert(fid);
    fprintf(fid,"%s", str);
  fclose(fid);