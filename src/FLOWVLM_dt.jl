# ------------ DATA TYPES ------------------------------------------------------
# Number types: Bugs can be easily cought during development by specifying all
# input types, however, in order to get automatic gradients types must be
# dynamic. The following flag turns number inputs into hard types (floats or
# ints) if true, or into dynamic types if false.
const dev_flag = false

if dev_flag
  const FWrap = Float64             # Float wrapper
  const IWrap = Int64               # Int wrapper
  const FArrWrap = Array{FWrap,1}   # Float array wrapper
  const IArrWrap = Array{IWrap,1}   # Int array wrapper
  const FMWrap = Array{FWrap,2}     # Float matrix wrapper
  const IMWrap = Array{IWrap,2}     # Int matrix wrapper
else
  const FWrap = Real
  const IWrap = Int64
  const FArrWrap = Array{T,1} where {T<:FWrap}
  const IArrWrap = Array{IWrap,1}
  const FMWrap = Array{T,2} where {T<:FWrap}
  const IMWrap = Array{IWrap,2}
end
