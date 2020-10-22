class RW_FOR:

  def wr_for(self,arr_in,fmt='%13.6E',n_lin=6):

    arr_flat=arr_in.transpose().ravel()
    nx=len(arr_flat)
    out_str=''
    for jx in range(nx):
      out_str +=  (fmt %arr_flat[jx])
      if (jx%n_lin == n_lin-1):
        out_str += '\n'
# Line break after writing data, but no double break
    if (nx%n_lin != 0):
      out_str += '\n'

    return out_str

  def ssplit(self,ll):

    nE = ll.count('E')
    lE = ll.rsplit('E')
    for i in range(nE):
      dd = lE[i]
      cc = lE[i+1]
      c = cc[0:3]
      if i == 0:
        d = dd
      else:
        d = dd[3:len(dd)]
      res = float(d+'E'+c)
      if i == 0:
        a = [res]
      else:
        a.append(res)

    return a
