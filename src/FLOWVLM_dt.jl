# ------------ DATA TYPES ------------------------------------------------------
# Number types: Bugs can be easily cought during development by specifying all
# input types, however, in order to get automatic gradients types must be
# dynamic. The following flag turns number inputs into hard types (floats or
# ints) if true, or into dynamic types if false.
const dev_flag = true

if dev_flag
  const FWrap = Float64             # Float wrapper
  const IWrap = Int64               # Int wrapper
  const IArrWrap = Array{IWrap,1}   # Int array wrapper
  const IMWrap = Array{IWrap,2}     # Int matrix wrapper
else
  const FWrap = Number
  const IWrap = Int64
  const IArrWrap = Array{IWrap,1}
  const IMWrap = AbstractArray{IWrap,2}
end
