function [] = write_tomo_multipleSection(folder_name,cube)

mkdir(folder_name)
for section = 1:size(cube,3)
   imwrite(uint8(cube(:,:,section)),[folder_name '/' num2str(section,'%05.f') '.tiff'],'TIFF') 
end