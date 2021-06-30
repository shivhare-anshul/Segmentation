function [Seg  Ncut] = Partition(I, W, segNcut, segArea)
N = length(W);
d = sum(W, 2);
D = spdiags(d, 0, N, N); 
[U,S] = eigs(D-W, D, 2, 'sm');
U2 = U(:, 2);
t = mean(U2);
t = fminsearch('Value', t, [], U2, W, D);
A = find(U2 > t);
B = find(U2 <= t);
ncut = Value(t, U2, W, D);
if (length(A) < segArea || length(B) < segArea) || ncut > segNcut
    Seg{1}   = I;
    Ncut{1} = ncut;       
    return;
end
[SegA  NcutA] = Partition(I(A), W(A, A), segNcut, segArea);
[SegB  NcutB] = Partition(I(B), W(B, B), segNcut, segArea);
Seg   = [SegA SegB];
Ncut = [NcutA NcutB];
end
