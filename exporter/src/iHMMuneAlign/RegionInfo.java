package iHMMuneAlign;


public class RegionInfo
{

    int start_pos;
    int end_pos;
    String region_string;

    public RegionInfo(String regionString, int startPos, int endPos)
    {
        region_string = regionString;
        start_pos = startPos;
        end_pos = endPos;
    }
}
