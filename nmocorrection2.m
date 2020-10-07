function seisnmo = nmocorrection2(t,dt,offset,seisdata,vnmo)
%nmocorrection2   NMO correction of CMP gather
%
%   seisnmo = nmo_correction2(t, dt, offset, seisdata, vnmo) applies NMO
%   correction to a CMP gather according the provided NMO velocities.
%
%   Input:
%   -t: Vector containing the travel times of the seismic data
%       (in seconds).
%   -dt: Number containing the time sampling of the seismic data
%        (in seconds).
%   -offset: Vector containing the source-receiver distance for each
%            seismic trace (in meters).
%   -seisdata: Matrix containing columns of seismic traces.
%   -vnmo: Vector containing NMO velocities (in meters/second).
%
%   Output:
%   -seisnmo: Matrix containing columns of NMO corrected seismic traces.

  seisnmo = zeros(size(seisdata));
  reflecttimes = sqrt(t.^2 + offset.^2./vnmo.^2);
  reflecttimes = reflecttimes.*(reflecttimes >= dt & reflecttimes <= t(end)-dt);
  for i=1:length(offset)
    seisnmo(:,i) = interp1(t,seisdata(:,i),reflecttimes(:,i),'spline',0);
  end
  
end