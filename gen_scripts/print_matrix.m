function str = print_matrix(X, type_str)
  if nargin < 2
    type_str = "double";
  end

  if type_str == "float"
    fmt_str = "%+1.9Ef";
    fmt_str_complex = "%+1.9Ef %+1.9Efi";
  elseif type_str == "double"
    fmt_str = "%+1.19E";
    fmt_str_complex = "%+1.19E %+1.19Ei";
  else
    error("type_str is unknown [float, double]")
  end

  str = "";
  if size(X,1) == 1
    str = str + sprintf("{ ");
  else
    str = str + sprintf("{\n    ");
  end
  for r_ = 1:size(X,1)
    for c_ = 1:size(X,2)
      if (abs(imag(X(r_,c_))) > 1E-9)
        str = str + sprintf(fmt_str_complex, real(X(r_,c_)), imag(X(r_,c_)));
      else
        str = str + sprintf(fmt_str, X(r_,c_));
      end
      if c_ < size(X,2)
        str = str + sprintf(", ");
      end
    end
    if r_ <  size(X,1)
      str = str + sprintf(", //\n    ");
    end
  end
  if size(X,1) > 1
    str = str + sprintf(" //\n    ");
  end
  str = str + sprintf("}");
end
