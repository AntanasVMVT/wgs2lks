// See https://aka.ms/new-console-template for more information



var lat = 56.1;
var lon = 23.8;

Console.WriteLine($"orig:\t{lat} \t{lon}");

var m = GeoConverter.GeoToGrid(lat, lon, 2);
var n = GeoConverter.GridToGeo(m.east, m.north);
Console.WriteLine($"p1:\t{n.lat}  \t{n.lon}");
Console.WriteLine($"p2:\t{m.east}  \t{m.north}");


Console.WriteLine("");




