function plot_imag(y)
figure; 
subplot(1,2,1); surf(real(y),'FaceAlpha',0.5); title('Real part');
subplot(1,2,2); surf(imag(y),'FaceAlpha',0.5); title('Imag part');
end